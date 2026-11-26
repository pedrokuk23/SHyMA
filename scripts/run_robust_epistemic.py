"""Settings of Concept D under the epistemic alternatives.

The selection under uncertainty (run_robust_selection.py) uses the base
dispersion model. Here the nominal setting of D and the settings selected by
the declared criteria are re-flown, with common random numbers (the
confirmation draws of run_robust_selection.py), under the epistemic
alternatives of the Monte Carlo analysis: booster calibration sampled from the
V01 propagation, from the leave-one-out spread and from the multi-flight
calibration, doubled engineering dispersions and the December NRLMSISE-00
atmosphere of CLA. Outputs: P(tau_use >= 3 s with the CLA link), its 95 %
interval and the 5th percentile of tau_use, and the paired difference of the
5th percentile between the selected and the nominal setting (bootstrap over
the common draws).
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
import run_robust_selection as rs  # noqa: E402
from shyma import dispersion as D  # noqa: E402

N = int(os.environ.get("RE_N", "500"))
SEED = rs.SEED_CONFIRM


def one(args):
    x, d = args
    kw = {"coast": x["coast"] + d["dcoast"]}
    if "isp" in d:
        kw.update(isp=d["isp"], s31_burn=d["b1"])
    atm = mc.msis_cla() if d.get("msis") else None
    v, tr = cs.fly("D", 400.0 * d["pl"], x["el"] + d["del"], az=cs.AZIMUTH + d["daz"], rho_factor=d["rho"],
                   dt=0.02, rec=5, wind=D.wind_of(d), dT=D.temp_of(d), beta_tv=cs.BETA_TV * d["beta"], beta_stack=cs.BETA_TV,
                   cd_scale=d["cd"], impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                   burn_factor=(d["burn1"], d["burn2"]), atm=atm, **kw)
    w = cs.windows("D", v, tr, cs.DSEP)[0]
    return cs.useful_window(tr, w, ant_loss_db=mc.ANT_LOSS_DB, fade_db=mc.FADE_DB, mask_deg=mc.ELEV_MASK_DEG)


def draws(kind, k):
    rng = random.Random(SEED)
    rc = random.Random(SEED + 2000)
    ds = [D.sample(rng, k, 2.0 if kind == "doubled_sigma" else 1.0) for _ in range(N)]
    cal = {"calibration": mc.calibration_uncertainty, "calibration_loo": mc.cal_loo,
           "calibration_multiflight": mc.cal_multiflight}.get(kind)
    if cal is not None:
        c = cal()
        for d in ds:
            z1, z2 = rc.gauss(0, 1), rc.gauss(0, 1)
            d["isp"] = c["isp"] + c["sigma_isp"] * z1
            d["b1"] = min(max(c["b1"] + c["sigma_b1"] * (c["corr"] * z1 + math.sqrt(1 - c["corr"] ** 2) * z2), 4.0), 20.0)
            d["k_wind"] = max(c["k_wind"] + c["sigma_k_wind"] * rc.gauss(0, 1), 0.0)
    if kind in ("msis", "calibration_multiflight"):   # the ten-flight calibration was fitted with NRLMSISE-00
        for d in ds:
            d["msis"] = True
    return ds


def p05_diff_ci(a, b, n_boot=1000, seed=3):
    r = random.Random(seed)
    n = len(a)
    vals = []
    for _ in range(n_boot):
        idx = [r.randrange(n) for _ in range(n)]
        vals.append(D.percentile([a[i] for i in idx], 0.05) - D.percentile([b[i] for i in idx], 0.05))
    return [D.percentile(vals, 0.025), D.percentile(vals, 0.975)]


if __name__ == "__main__":
    t0 = time.time()
    sel = json.load(open(os.path.join(rs.ROOT, "results", "robust_selection.json")))["concepts"]["D"]
    settings = {"nominal": sel["nominal_x"], "R1": sel["R1_x"], "R3": sel["R3_x"], "primary": sel["primary_x"]}
    uniq = []
    for x in settings.values():
        if x and x not in uniq:
            uniq.append(x)
    k = mc.k_wind_calibrated()
    out = {"n": N, "seed": SEED, "settings": settings, "sets": {}}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as pool:
        kinds = ("base", "calibration", "calibration_loo", "calibration_multiflight", "doubled_sigma", "msis")
        only = [k for k in os.environ.get("EPI_SETS", "").split(",") if k]
        if only:   # recompute some sets only, keeping the others of the existing result file
            prev = json.load(open(os.path.join(rs.ROOT, "results", "robust_epistemic.json")))
            out["sets"].update({k: v for k, v in prev["sets"].items() if k not in only})
        for kind in (k for k in kinds if not only or k in only):
            ds = draws(kind, k)
            res = {}
            taus = {}
            for x in uniq:
                tu = pool.map(one, [(x, d) for d in ds], chunksize=4)
                kk = sum(t >= rs.TAU_REQ for t in tu)
                key = "%.0f/%.0f" % (x["el"], x["coast"])
                taus[key] = tu
                res[key] = {"x": x, "p": kk / N, "p_ci95": D.clopper_pearson(kk, N),
                            "tau_use_p05": D.percentile(tu, 0.05), "tau_use_p50": D.percentile(tu, 0.5)}
                print(kind, key, round(kk / N, 3), round(res[key]["tau_use_p05"], 2), flush=True)
            nk = "%.0f/%.0f" % (settings["nominal"]["el"], settings["nominal"]["coast"])
            pk = "%.0f/%.0f" % (settings["primary"]["el"], settings["primary"]["coast"])
            res["primary_minus_nominal_p05"] = D.percentile(taus[pk], 0.05) - D.percentile(taus[nk], 0.05)
            res["primary_minus_nominal_p05_ci95"] = p05_diff_ci(taus[pk], taus[nk])
            out["sets"][kind] = res
    out["cpu_s"] = time.time() - t0
    json.dump(out, open(os.path.join(rs.ROOT, "results", "robust_epistemic.json"), "w"), indent=1)
