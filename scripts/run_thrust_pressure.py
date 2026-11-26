"""Sensitivity of the case study to the ambient-pressure thrust loss.

The baseline model uses F = mdot Isp_vac g0 - p_amb A_e with an exit diameter
of 0.45 m for both motors (vehicle.EXIT_DIAMETER_M; no public nozzle data), with
Isp_vac and the S31 burn time calibrated on flight V01. Here the exit diameter
is varied (0: constant equivalent thrust, the model of earlier versions; 0.40
and 0.50 m), the booster is recalibrated on V01 for each variant, and the
nominal concepts are flown again at their baseline settings.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
import run_vsb30_validation as vv  # noqa: E402
from scipy.stats import chi2  # noqa: E402
from shyma import calibration, vehicle  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
EXIT_DIAMETERS_M = (0.0, 0.40, 0.50)


def held_out(isp, b1):
    """Held-out VSB-30 apogee statistics with the current thrust model: RMS and mean error (%),
    and the upper-tail chi-square test of the residuals normalised by the Monte Carlo flight
    dispersion of results/vsb30_validation.json."""
    vj = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    sg = {f["flight"]: f["mc_sigma_apogee_km"] for f in vj["flights"]}
    e, z = [], []
    for r in vv.load():
        if r["flight"] == "V01":
            continue
        a = vv.predict(r, isp, b1)["apogee_km"]
        act = float(r["apogee_actual_km"])
        e.append(100 * (a - act) / act)
        z.append((a - act) / sg[r["flight"]])
    q = sum(x * x for x in z)
    return {"apogee_rms_pct": math.sqrt(sum(x * x for x in e) / len(e)), "apogee_mean_pct": sum(e) / len(e),
            "chi2": q, "p_upper": float(chi2.sf(q, len(z)))}


def nominal_set():
    out = {}
    for c in "ABCD":
        v, tr = cs.fly(c, 400.0, cs.NOMINAL[c])
        s = cs.summarize(c, v, tr, light=True)
        out[c] = {k: s[k] for k in ("window_s", "burnout_mach", "apogee_km", "impact_range_km")}
    return out


if __name__ == "__main__":
    t0 = time.time()
    isp0, b10 = vehicle.ISP_VSB30, vehicle.S31_BURN
    ae0 = vehicle.EXIT_AREA
    res = {"exit_diameters_m": EXIT_DIAMETERS_M, "baseline": nominal_set(), "baseline_isp_s": isp0,
           "baseline_held_out": held_out(isp0, b10),
           "baseline_exit_diameter_m": vehicle.EXIT_DIAMETER_M,
           "baseline_s31_burn_s": b10, "cases": {}}
    for de in EXIT_DIAMETERS_M:
        ae = math.pi * de ** 2 / 4.0
        vehicle.EXIT_AREA = (ae, ae)
        c = calibration.calibrate(calibration.V01, x0=(isp0 + 8.0 * (de / 0.45) ** 2 - 8.0, b10))
        vehicle.ISP_VSB30, vehicle.S31_BURN = c["isp"], c["b1"]
        res["cases"]["%.2f" % de] = {"exit_area_m2": ae, "isp_vac_s": c["isp"], "s31_burn_s": c["b1"],
                                     "concepts": nominal_set(), "held_out": held_out(c["isp"], c["b1"])}
        vehicle.ISP_VSB30, vehicle.S31_BURN = isp0, b10
        vehicle.EXIT_AREA = ae0
        print(de, res["cases"]["%.2f" % de], flush=True)
    res["cpu_s"] = time.time() - t0
    json.dump(res, open(os.path.join(ROOT, "results", "thrust_pressure.json"), "w"), indent=1)
