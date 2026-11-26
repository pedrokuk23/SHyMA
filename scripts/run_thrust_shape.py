"""Sensitivity of the Concept A and D dwell to the thrust-curve shape.

The paper's booster uses a constant equivalent thrust per motor. Here each
motor may instead follow a linear thrust law F = F_mean (1 + k (1 - 2 tau/t_b))
with the same impulse and burn time (k > 0 regressive, k < 0 progressive).
For every variant the booster is recalibrated on flight V01 (Isp and the
equivalent S31 burn time), so that all variants reproduce the calibration
apogee and range exactly; then Concepts A and D at their nominal settings
are flown at their nominal settings and with their launch settings re-tuned
(A: elevation; D: elevation and coast).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
from shyma import calibration, vehicle  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
K = 0.3
VARIANTS = {
    "constant": (0.0, 0.0),
    "S30 regressive": (0.0, K), "S30 progressive": (0.0, -K),
    "S31 regressive": (K, 0.0), "S31 progressive": (-K, 0.0),
    "both regressive": (K, K), "both progressive": (-K, -K),
}
DT = 0.02


def dwell(c, el, isp, b1, **kw):
    v, tr = cs.fly(c, 400.0, el, dt=DT, rec=1, isp=isp, s31_burn=b1, **kw)
    w = cs.windows(c, v, tr)[0]
    s30 = [s for s in v.stages if s.name == "S30"][0]
    k = min(range(len(tr.t)), key=lambda i: abs(tr.t[i] - v.t_burnout))
    return {"dwell_s": w["duration"], "burnout_mach": tr.mach[k], "burnout_alt_km": tr.alt[k] / 1e3,
            "s30_ign_s": s30.t_ign}


def main():
    t0 = time.time()
    out = {"k": K, "dt": DT, "variants": {}}
    for name, shape in VARIANTS.items():
        vehicle.THRUST_SHAPE = shape
        cal = calibration.calibrate(calibration.V01, x0=(vehicle.ISP_VSB30, vehicle.S31_BURN))
        isp, b1 = cal["isp"], cal["b1"]
        rec = {"shape": shape, "isp": isp, "b1": b1, "apogee_km": cal["apogee_km"], "range_km": cal["range_km"]}
        rec["A_nominal"] = dwell("A", cs.NOMINAL["A"], isp, b1)
        rec["D_nominal"] = dwell("D", cs.NOMINAL["D"], isp, b1, coast=cs.COAST_D)
        best_a = max(((el, dwell("A", el, isp, b1)["dwell_s"]) for el in [56 + 0.5 * i for i in range(25)]),
                     key=lambda x: x[1])
        rec["A_retuned"] = {"el": best_a[0], "dwell_s": best_a[1]}
        best_d = max(((el, co, dwell("D", el, isp, b1, coast=co)["dwell_s"])
                      for el in (69.0, 71.0, 73.0, 75.0, 77.0) for co in [17 + i for i in range(14)]),
                     key=lambda x: x[2])
        rec["D_retuned"] = {"el": best_d[0], "coast": best_d[1], "dwell_s": best_d[2]}
        out["variants"][name] = rec
        print(name, json.dumps({k: rec[k] for k in ("isp", "b1", "A_nominal", "D_nominal", "A_retuned", "D_retuned")}),
              flush=True)
    vehicle.THRUST_SHAPE = (0.0, 0.0)
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "thrust_shape.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
