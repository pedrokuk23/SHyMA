"""Residual dispersion represented as a launch tip-off instead of a residual wind.

The lumped residual-dispersion parameter k_w acts through the wind, and at the
depressed elevations of Concepts A and D its effect differs from that of the
sources it represents (thrust misalignment, tip-off), which act on the attitude.
Here an equivalent tip-off attitude error just after rail exit (two-axis
Gaussian, standard deviation sigma_t per axis, uniform phase) replaces k_w:
sigma_t is calibrated on V01 with k_w = 0 to the same 3-sigma radius of 50 km,
and the case-study Monte Carlo of A, B and D is repeated with it.
"""
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
import run_montecarlo as mc  # noqa: E402
from shyma import dispersion as D  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
N_CAL = int(os.environ.get("MC_N_TIP_CAL", "400"))
N_CS = int(os.environ.get("MC_N_TIP", "500"))
SIG_GRID = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
CFG = (392.0, 84.0, 45.0, "CLA")
TARGET_R = 50.0
TAU_REQ = 3.0


def add_tip(d, rng, sig):
    a, b = rng.gauss(0, 1), rng.gauss(0, 1)
    d = dict(d, k_wind=0.0)
    if sig > 0:
        d["tip"] = (sig * math.hypot(a, b), math.atan2(b, a))
    return d


def r3(out):
    n = len(out)
    lats, lons = [x["impact_lat"] for x in out], [x["impact_lon"] for x in out]
    return D.ellipse(sum(lats) / n, sum(lons) / n, lats, lons)["r3sig_gauss_km"]


def one(args):
    c, d = args
    kw = {"coast": cs.COAST_D + d["dcoast"]} if c == "D" else {}
    if "tip" in d:
        kw["tip_att"] = d["tip"]
    v, tr = cs.fly(c, 400.0 * d["pl"], cs.NOMINAL[c] + d["del"], az=cs.AZIMUTH + d["daz"], rho_factor=d["rho"],
                   dt=0.02, rec=5, wind=None, dT=D.temp_of(d), beta_tv=cs.BETA_TV * d["beta"], beta_stack=cs.BETA_TV,
                   cd_scale=d["cd"], impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                   burn_factor=(d["burn1"], d["burn2"]), **kw)
    w = cs.windows(c, v, tr, max(cs.DSEP + d["dsep"], 0.0))[0]
    return {"met": bool(w["met"]), "window_s": w["duration"],
            "tau_use": cs.useful_window(tr, w, ant_loss_db=mc.ANT_LOSS_DB, fade_db=mc.FADE_DB,
                                        mask_deg=mc.ELEV_MASK_DEG),
            "impact_lat": tr.lat[-1], "impact_lon": tr.lon[-1]}


if __name__ == "__main__":
    t0 = time.time()
    rng = random.Random(21)
    base = [D.sample(rng, 0.0) for _ in range(N_CAL)]
    rt = random.Random(22)
    unit = [(rt.gauss(0, 1), rt.gauss(0, 1)) for _ in range(N_CAL)]
    curve = {}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as p:
        for s in SIG_GRID:
            dr = []
            for d, (a, b) in zip(base, unit):
                d = dict(d, k_wind=0.0)
                if s > 0:
                    d["tip"] = (s * math.hypot(a, b), math.atan2(b, a))
                dr.append((CFG, d))
            curve[s] = r3(p.map(D.one_flight, dr))
            print("sig", s, curve[s], flush=True)
        sig = float("nan")
        for a, b in zip(SIG_GRID[:-1], SIG_GRID[1:]):
            if (curve[a] - TARGET_R) * (curve[b] - TARGET_R) <= 0 and curve[a] != curve[b]:
                sig = a + (TARGET_R - curve[a]) * (b - a) / (curve[b] - curve[a])
                break
        if not sig == sig:
            raise RuntimeError("tip-off sigma not bracketed by SIG_GRID %s (3-sigma radii %s)" % (SIG_GRID, curve))
        res = {"sigma_tip_deg": sig, "curve": {str(k): v for k, v in curve.items()}, "N_cal": N_CAL,
               "N": N_CS, "concepts": {}}
        for c in "ABD":
            r = random.Random(31 + ord(c))
            draws = [add_tip(D.sample(r, 0.0), r, sig) for _ in range(N_CS)]
            out = p.map(one, [(c, d) for d in draws])
            k = sum(o["met"] for o in out)
            ku = sum(o["tau_use"] >= TAU_REQ for o in out)
            res["concepts"][c] = {"p_window": k / N_CS, "p_window_ci95": D.clopper_pearson(k, N_CS),
                                  "p_use": ku / N_CS, "p_use_ci95": D.clopper_pearson(ku, N_CS),
                                  "window_p05": D.percentile([o["window_s"] for o in out], 0.05),
                                  "tau_use_p05": D.percentile([o["tau_use"] for o in out], 0.05),
                                  "r3_km": r3(out)}
            print(c, res["concepts"][c], flush=True)
    res["wall_s"] = time.time() - t0
    json.dump(res, open(os.path.join(ROOT, "results", "tipoff_check.json"), "w"), indent=1)
