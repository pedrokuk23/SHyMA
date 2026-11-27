"""Selection of the launch settings of each concept under uncertainty.

The deterministic inverse design (run_optimizer.py) and the Monte Carlo
analysis of the nominal settings (run_montecarlo.py) are joined here: a grid
of candidate settings of each architecture (the nominal setting included) is
evaluated with the same dispersed draws (common random numbers), the
candidates are ranked by a declared probabilistic criterion, and the
selected setting is re-evaluated, together with the nominal one, on an
independent sample of the dispersions.

Measure of effectiveness: the useful test time tau_use, the longest
continuous interval in which the threshold window is met and the CLA link
closes (degraded margin >= 0 dB, elevation >= 3 deg; run_casestudy.
useful_window). Disjoint pieces are not added. The same quantity without the
link condition (on-board recording) is also kept.

Requirement: tau_use >= TAU_REQ (3 s, the minimum duration of the threshold
window). Criteria:
  R1  maximise the lower 95 % Clopper-Pearson bound of P(tau_use >= TAU_REQ);
      ties broken by the 5th percentile of tau_use;
  R2  the shortest mean impact range among the candidates whose lower bound
      is >= P_MIN;
  R3  the largest 5th percentile of tau_use (the setting that keeps the most
      useful time when conditions vary);
  primary (declared criterion of the paper): the largest 5th percentile of
      tau_use among the candidates whose lower bound is within PRIMARY_TOL of
      the best lower bound (probability first, then robustness of the time).
Probabilities are conditional on the dispersion model of shyma.dispersion
(model alternatives such as the thrust shape are scenarios, not samples).
The screening maximum over many candidates is optimistic (selection bias);
the confirmation sample removes it for the selected setting.
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
N_SCREEN = int(os.environ.get("RS_N", "200"))
N_CONFIRM = int(os.environ.get("RS_NC", "500"))
SEED_SCREEN = 20261023
SEED_CONFIRM = 20261024
TAU_REQ = 3.0
P_MIN = 0.90
PRIMARY_TOL = 0.01   # tolerance on the lower 95 % bound of P for the primary criterion
TAU_GRID = [0.5 * i for i in range(41)]

CANDIDATES = {
    "A": [{"el": e} for e in (57.0, 59.0, 61.0, 63.0, 65.0, 67.0, 69.0)],
    "B": [{"el": e} for e in (60.0, 66.0, 72.0, 78.0, 84.0, 85.0, 86.0, 87.0, 88.0)],
    "C": [{"el": e} for e in (60.0, 66.0, 72.0, 78.0, 84.0, 85.0, 86.0, 87.0, 88.0)],
    "D": [{"el": e, "coast": c} for e in (68.5, 70.5, 72.5, 74.5, 76.5) for c in (12.0, 15.0, 18.0, 21.0, 24.0, 27.0)],
}


def nominal(c):
    x = {"el": cs.NOMINAL[c]}
    if c == "D":
        x["coast"] = cs.COAST_D
    return x


def one(args):
    c, x, d = args
    kw = {"coast": x["coast"] + d["dcoast"]} if c == "D" else {}
    v, tr = cs.fly(c, 400.0 * d["pl"], x["el"] + d["del"], az=cs.AZIMUTH + d["daz"], rho_factor=d["rho"],
                   dt=0.02, rec=5, wind=D.wind_of(d), dT=D.temp_of(d), beta_tv=cs.BETA_TV * d["beta"], beta_stack=cs.BETA_TV,
                   cd_scale=d["cd"], impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                   burn_factor=(d["burn1"], d["burn2"]), **kw)
    dsep = max(cs.DSEP + d["dsep"], 0.0)
    w = cs.windows(c, v, tr, dsep)[0]
    return {"tau_window": w["duration"] if w["met"] or w["duration"] > 0 else 0.0,
            "tau_use": cs.useful_window(tr, w, ant_loss_db=mc.ANT_LOSS_DB, fade_db=mc.FADE_DB,
                                        mask_deg=mc.ELEV_MASK_DEG),
            "tau_rec": cs.useful_window(tr, w, need_link=False),
            "range_km": cs.analysis.downrange_km(tr, cs.SITE)[-1]}


def stats(out):
    n = len(out)
    tu = [o["tau_use"] for o in out]
    k = sum(t >= TAU_REQ for t in tu)
    kr = sum(o["tau_rec"] >= TAU_REQ for o in out)
    lo, hi = D.clopper_pearson(k, n)
    curve = []
    for tq in TAU_GRID:
        m = sum(t >= tq for t in tu)
        curve.append([tq, m / n, D.clopper_pearson(m, n)[0]])
    return {"n": n, "p_req": k / n, "p_req_ci95": [lo, hi], "p_req_rec": kr / n,
            "tau_use_p05": D.percentile(tu, 0.05), "tau_use_p50": D.percentile(tu, 0.5),
            "tau_window_p50": D.percentile([o["tau_window"] for o in out], 0.5),
            "range_mean_km": sum(o["range_km"] for o in out) / n, "curve": curve}


def draws(seed, n):
    rng = random.Random(seed)
    k = mc.k_wind_calibrated()
    return [D.sample(rng, k, 1.0) for _ in range(n)]


def evaluate(pool, c, x, ds):
    return stats(pool.map(one, [(c, x, d) for d in ds], chunksize=4))


def main():
    t0 = time.time()
    ds = draws(SEED_SCREEN, N_SCREEN)
    dc = draws(SEED_CONFIRM, N_CONFIRM)
    out = {"tau_req": TAU_REQ, "p_min": P_MIN, "n_screen": N_SCREEN, "n_confirm": N_CONFIRM,
           "ant_loss_db": mc.ANT_LOSS_DB, "fade_db": mc.FADE_DB, "concepts": {}}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as pool:
        for c, cands in CANDIDATES.items():
            nom = nominal(c)
            assert nom in cands
            rows = []
            for x in cands:
                s = evaluate(pool, c, x, ds)
                rows.append({"x": x, **s})
                print(c, x, "P=%.3f [%.3f] p05=%.2f p50=%.2f R=%.0f" % (
                    s["p_req"], s["p_req_ci95"][0], s["tau_use_p05"], s["tau_use_p50"], s["range_mean_km"]),
                    flush=True)
            r1 = max(rows, key=lambda r: (r["p_req_ci95"][0], r["tau_use_p05"]))
            ok = [r for r in rows if r["p_req_ci95"][0] >= P_MIN]
            r2 = min(ok, key=lambda r: r["range_mean_km"]) if ok else None
            env = [[tq, max(r["curve"][i][2] for r in rows)] for i, tq in enumerate(TAU_GRID)]
            r3 = max(rows, key=lambda r: (r["tau_use_p05"], r["p_req_ci95"][0]))
            # primary criterion (declared in the paper): the largest 5th percentile of tau_use among the
            # candidates whose lower 95 % bound of P is within PRIMARY_TOL of the best lower bound
            best_lcb = max(r["p_req_ci95"][0] for r in rows)
            rp = max([r for r in rows if r["p_req_ci95"][0] >= best_lcb - PRIMARY_TOL],
                     key=lambda r: (r["tau_use_p05"], r["p_req_ci95"][0]))
            rec = {"screen": rows, "nominal_x": nom, "R1_x": r1["x"], "R2_x": r2["x"] if r2 else None,
                   "R3_x": r3["x"], "primary_x": rp["x"], "primary_tol": PRIMARY_TOL,
                   "envelope_lcb": env, "confirm": {}}
            for tag, x in (("nominal", nom), ("R1", r1["x"]), ("R2", r2["x"] if r2 else None), ("R3", r3["x"]),
                           ("primary", rp["x"])):
                if x is None:
                    continue
                same = [t for t, v in rec["confirm"].items() if v["x"] == x]
                if same:
                    rec["confirm"][tag] = rec["confirm"][same[0]]
                    continue
                s = evaluate(pool, c, x, dc)
                rec["confirm"][tag] = {"x": x, **s}
                print("confirm", c, tag, x, "P=%.3f %s p05=%.2f p50=%.2f R=%.0f" % (
                    s["p_req"], s["p_req_ci95"], s["tau_use_p05"], s["tau_use_p50"], s["range_mean_km"]),
                    flush=True)
            out["concepts"][c] = rec
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "robust_selection.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
