"""Local screen around the setting of Concept D selected under uncertainty.

The grid surrounds the setting selected in run_robust_selection.py by the
primary criterion: elevation from 4 deg below to 2 deg above in 1 deg steps
and coast from 6 s below to 3 s above in 1.5 s steps. All points use the same
screening draws and requirement (tau_use >= 3 s with a closed CLA link). To
avoid the winner's bias of comparing the best screened point on the draws that
selected it, the best grid point and the selected setting are then re-flown on
the independent confirmation draws, and the difference of their 5th
percentiles is given with a paired bootstrap interval.
"""
import json
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(__file__))
import run_robust_selection as rs  # noqa: E402
from shyma import dispersion as D  # noqa: E402

N_BOOT = 1000


def grid(sel):
    e0, c0 = sel["el"], sel["coast"]
    return [{"el": e0 + de, "coast": max(c0 + dc, 0.0)} for de in (-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0)
            for dc in (-6.0, -4.5, -3.0, -1.5, 0.0, 1.5, 3.0)]


if __name__ == "__main__":
    t0 = time.time()
    sel = json.load(open(os.path.join(rs.ROOT, "results", "robust_selection.json")))["concepts"]["D"]["primary_x"]
    ds = rs.draws(rs.SEED_SCREEN, rs.N_SCREEN)
    out = {"n_screen": rs.N_SCREEN, "seed": rs.SEED_SCREEN, "tau_req": rs.TAU_REQ, "selected": sel, "grid": []}
    taus = {}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as pool:
        for x in grid(sel):
            raw = pool.map(rs.one, [("D", x, d) for d in ds], chunksize=4)
            st = rs.stats(raw)
            st.pop("curve", None)
            out["grid"].append({"x": x, **st})
            taus[(x["el"], x["coast"])] = [o["tau_use"] for o in raw]
            print(x, round(st["p_req"], 3), round(st["tau_use_p05"], 2), round(st["tau_use_p50"], 2),
                  round(st["range_mean_km"]), flush=True)
        best = max(out["grid"], key=lambda r: (r["tau_use_p05"], r["p_req_ci95"][0]))
        dc = rs.draws(rs.SEED_CONFIRM, rs.N_CONFIRM)
        tb = [o["tau_use"] for o in pool.map(rs.one, [("D", best["x"], d) for d in dc], chunksize=4)]
        ts = [o["tau_use"] for o in pool.map(rs.one, [("D", sel, d) for d in dc], chunksize=4)]
    r = random.Random(5)
    n = len(ts)
    diffs = []
    for _ in range(N_BOOT):
        idx = [r.randrange(n) for _ in range(n)]
        diffs.append(D.percentile([tb[i] for i in idx], 0.05) - D.percentile([ts[i] for i in idx], 0.05))
    out["best"] = best["x"]
    out["confirm_n"] = n
    out["best_confirm_p05"], out["selected_confirm_p05"] = D.percentile(tb, 0.05), D.percentile(ts, 0.05)
    out["best_confirm_p"] = sum(t >= rs.TAU_REQ for t in tb) / n
    out["selected_confirm_p"] = sum(t >= rs.TAU_REQ for t in ts) / n
    out["best_minus_selected_p05"] = out["best_confirm_p05"] - out["selected_confirm_p05"]
    out["best_minus_selected_p05_ci95"] = [D.percentile(diffs, 0.025), D.percentile(diffs, 0.975)]
    es, cs_ = [g["x"]["el"] for g in out["grid"]], [g["x"]["coast"] for g in out["grid"]]
    out["best_on_edge"] = best["x"]["el"] in (min(es), max(es)) or best["x"]["coast"] in (min(cs_), max(cs_))
    out["cpu_s"] = time.time() - t0
    print(out["best"], out["best_minus_selected_p05"], out["best_minus_selected_p05_ci95"], out["best_on_edge"])
    json.dump(out, open(os.path.join(rs.ROOT, "results", "robust_local.json"), "w"), indent=1)
