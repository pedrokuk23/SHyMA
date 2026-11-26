"""Geometry-based zero-lift drag of the VSB-30 (shyma.drag) and its effect on
the mission results, alone and together with the geometry-based stability
derivatives (shyma.stability), compared with the generic drag curve of the paper.

Variants: geometric drag (3 fins, t/c = 0.03, 3 mm leading-edge radius, nozzle
exit / base area 0.72 from the scale model); thicker fins (t/c = 0.05); smaller
nozzle exit (area ratio 0.5, larger powered base drag); and geometric drag with
geometric stability derivatives. For each variant the booster is recalibrated on
V01, then the held-out VSB-30 apogee errors, the unit-wind effect and the dwell
of Concepts A and D are computed (run_aero_sensitivity.py functions).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import drag as D, flight, stability as St, vehicle  # noqa: E402
import run_casestudy as cs  # noqa: E402
import run_geometric_aero as rg  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
MACHS = [0.3, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
RE_PER_M = 5e6


def retune(isp, b1, sfac):
    """Launch settings re-tuned for the variant (grids of run_thrust_shape.py)."""
    flight.AERO_FACTORY = sfac
    try:
        def dw(c, el, **kw):
            v, tr = cs.fly(c, 400.0, el, dt=0.02, rec=1, isp=isp, s31_burn=b1, **kw)
            return cs.windows(c, v, tr)[0]["duration"]
        ra_ = max(((el, dw("A", el)) for el in [56 + 0.5 * i for i in range(25)]), key=lambda z: z[1])
        rd = max(((el, co, dw("D", el, coast=co)) for el in (69.0, 71.0, 73.0, 75.0, 77.0)
                  for co in [14 + i for i in range(15)]), key=lambda z: z[2])
    finally:
        flight.AERO_FACTORY = None
    return {"A_retuned": {"el": ra_[0], "dwell_s": ra_[1]}, "D_retuned": {"el": rd[0], "coast": rd[1], "dwell_s": rd[2]}}


def table(model):
    rows = []
    for M in MACHS:
        full_off = model.components(0, M, RE_PER_M, False)
        full_on = model.components(0, M, RE_PER_M, True)
        up_off = model.components(1, M, RE_PER_M, False)
        up_on = model.components(1, M, RE_PER_M, True)
        rows.append({"mach": M, "full_off": full_off, "full_on": full_on["total"], "upper_off": up_off["total"],
                     "upper_on": up_on["total"], "generic": vehicle.cd_generic(M)})
    return rows


if __name__ == "__main__":
    t0 = time.time()
    base = D.vsb30_drag()
    out = {"re_per_m": RE_PER_M, "table": table(base), "variants": {}}
    variants = {
        "geometric drag": (lambda v: D.attach(v, base), None),
        "geometric drag, t/c 0.05": (lambda v: D.attach(v, D.vsb30_drag(t_c=0.05)), None),
        "geometric drag, exit ratio 0.5": (lambda v: D.attach(v, D.vsb30_drag(exit_area_ratio=0.5)), None),
        "geometric drag and stability": (lambda v: D.attach(v, base), lambda v: St.vsb30_aero(v, 3)),
    }
    for name, (dfac, sfac) in variants.items():
        vehicle.DRAG_FACTORY = dfac
        try:
            out["variants"][name] = rg.mission(sfac)
            out["variants"][name].update(retune(out["variants"][name]["isp"], out["variants"][name]["b1"], sfac))
        finally:
            vehicle.DRAG_FACTORY = None
        print(name, json.dumps(out["variants"][name]), flush=True)
    for cc in "AD":
        d = [r["dwell_" + cc] for r in out["variants"].values()]
        out["dwell_range_" + cc] = [min(d), max(d)]
        d = [r[cc + "_retuned"]["dwell_s"] for r in out["variants"].values()]
        out["dwell_retuned_range_" + cc] = [min(d), max(d)]
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "geometric_drag.json"), "w") as f:
        json.dump(out, f, indent=1)
