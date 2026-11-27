"""Sensitivity of the results to the VSB-30 stage masses.

The stage masses of the model (vehicle.S31_PROP/S31_DRY, S30_PROP/S30_DRY) combine
the IAE motor table with an assumed split of the inert parts. This script scales the
dry masses (+-10 %) and the propellant masses (+-5 %) of both stages,
recalibrates the specific impulse and the S31 burn time on flight V01 for each
case, and reports the held-out apogee error on the other nine VSB-30 flights
and the nominal threshold-window dwell and impact range of the four test
concepts. Output: results/stage_mass_sensitivity.json.
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from shyma import calibration, vehicle  # noqa: E402
import run_vsb30_validation as rv  # noqa: E402
import run_casestudy as rc  # noqa: E402

CASES = [("base", 1.0, 1.0), ("dry -10%", 1.0, 0.9), ("dry +10%", 1.0, 1.1),
         ("propellant -5%", 0.95, 1.0), ("propellant +5%", 1.05, 1.0)]
ORIG = vehicle.vsb30
ISP0, B10 = vehicle.ISP_VSB30, vehicle.S31_BURN


def patched(pf, df):
    def f(*a, **k):
        k.setdefault("impulse_factor", (pf, pf))
        k.setdefault("dry_factor", df)
        return ORIG(*a, **k)
    return f


def main():
    t0 = time.time()
    rows = rv.load()
    out = {"cases": [], "concepts": "ABCD", "note": "Isp and t_b1 recalibrated on V01 for each case"}
    for label, pf, df in CASES:
        vehicle.vsb30 = patched(pf, df)
        c = calibration.calibrate()
        vehicle.ISP_VSB30, vehicle.S31_BURN = c["isp"], c["b1"]
        recs = rv.held_out_set(rows, c["isp"], c["b1"])
        s = rv.summary(recs)["apogee_model_pct"]
        conc = {}
        for k in "ABCD":
            v, tr = rc.fly(k, 400.0, rc.NOMINAL[k], rec=1)
            sm = rc.summarize(k, v, tr)
            conc[k] = {"window_s": sm["window_s"], "impact_range_km": sm.get("impact_range_km")}
        out["cases"].append({"case": label, "prop_factor": pf, "dry_factor": df,
                             "isp": c["isp"], "s31_burn": c["b1"],
                             "heldout_apogee_rms_pct": s["rms"], "heldout_apogee_mean_pct": s["mean"],
                             "concepts": conc})
        print(label, round(c["isp"], 2), round(c["b1"], 3), round(s["rms"], 2),
              {k: round(conc[k]["window_s"], 2) for k in conc}, flush=True)
    vehicle.vsb30, vehicle.ISP_VSB30, vehicle.S31_BURN = ORIG, ISP0, B10
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "stage_mass_sensitivity.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
