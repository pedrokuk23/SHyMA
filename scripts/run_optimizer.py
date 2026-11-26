"""Inverse design of the launch settings (scripts/launch_optimizer.py).

1) For each concept and a 400 kg test vehicle, the launch elevation (and the
   interstage coast for Concept D) that maximises the threshold dwell.
2) For a 350 kg vehicle, the settings that maximise the objective dwell.
3) Mass ceiling: the largest test-vehicle mass (10 kg resolution) for which
   some setting meets the threshold window, found by bisection on the mass
   with an optimisation at every step; the closest setting and the violated
   requirements are reported just above the ceiling.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
import launch_optimizer as lo  # noqa: E402
import run_casestudy as cs     # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
FREE = {"A": {"el": (54, 88, 2)}, "B": {"el": (54, 88, 2)}, "C": {"el": (54, 88, 2)},
        "D": {"el": (54, 88, 2), "coast": (0, 30, 3)}}


def fixed(c, m):
    return {"payload": m, "az": cs.AZIMUTH, "beta": cs.BETA_TV, "el": cs.NOMINAL[c], "coast": cs.COAST_D,
            "dsep": cs.DSEP}


def brief(r):
    b = r["best"]
    return {"x": b["x"], "threshold_s": b["threshold"], "objective_s": b["objective"], "steady_s": b["steady"],
            "threshold_met": b["threshold_met"], "objective_met": b["objective_met"],
            "impact_range_km": b["impact_range_km"], "q_max_ascent_kpa": b["q_max_ascent_kpa"],
            "feasible": r["feasible"], "any_feasible": r["any_feasible"], "n_eval": r["n_eval"]}


def closest_report(c, m):
    r = lo.search(c, fixed(c, m), FREE[c], "threshold", "max_dwell")
    x = dict(fixed(c, m), **r["best"]["x"])
    v, tr = cs.fly(c, m, x["el"], az=x["az"], beta_tv=x["beta"], dt=0.05, rec=2,
                   **({"coast": x["coast"]} if c == "D" else {}))
    w = lo.window_objects(c, v, x["dsep"])["threshold"]
    d = lo.diagnose(w, tr)
    return r, {"x": r["best"]["x"], "closest": d.get("closest"),
               "violated": [i["name"] for i in d["items"] if not i["ok"] and i["name"] != "Dwell"]}


if __name__ == "__main__":
    t0 = time.time()
    out = {"free": FREE, "dt": 0.05, "max_dwell_400": {}, "objective_350": {}, "mass_ceiling": {}}
    for c in "ABCD":
        out["max_dwell_400"][c] = brief(lo.search(c, fixed(c, 400.0), FREE[c], "threshold", "max_dwell"))
        # the optimum re-evaluated with the step of the case study (0.01 s), so that it
        # compares with the nominal dwell at the same numerical accuracy
        xf = dict(fixed(c, 400.0), **out["max_dwell_400"][c]["x"])
        rf = lo.quick(c, 400.0, xf["el"], xf["az"], xf["beta"], xf["coast"], xf["dsep"], "threshold", dt=0.01)
        out["max_dwell_400"][c]["threshold_s_fine"] = rf["threshold"]
        print("400", c, out["max_dwell_400"][c], flush=True)
        out["objective_350"][c] = brief(lo.search(c, fixed(c, 350.0), FREE[c], "objective", "max_dwell"))
        print("350", c, out["objective_350"][c], flush=True)
    for c in "ABCD":
        lo_m, hi_m = 400.0, 700.0
        r_hi, rep_hi = closest_report(c, hi_m)
        if r_hi["any_feasible"]:
            out["mass_ceiling"][c] = {"ceiling_kg": None, "note": "feasible at 700 kg"}
            continue
        while hi_m - lo_m > 10.0:
            mid = round((lo_m + hi_m) / 20.0) * 10.0
            r, rep = closest_report(c, mid)
            if r["any_feasible"]:
                lo_m = mid
            else:
                hi_m, rep_hi = mid, rep
            print("ceiling", c, lo_m, hi_m, flush=True)
        out["mass_ceiling"][c] = {"ceiling_kg": lo_m, "infeasible_kg": hi_m, "closest_at_infeasible": rep_hi}
        print("ceiling", c, out["mass_ceiling"][c], flush=True)
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "optimizer.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("wall", out["cpu_s"])
