"""Geometry-based stability derivatives of the VSB-30 (shyma.stability) and
their effect on the mission results, compared with the generic constants of
the paper (C_N_alpha = 12 /rad, static margin 1.5 calibers, C_mq = -300).

Geometry: shyma.geometry.vsb30_geometry (public dimensions and a 1:96 scale
model; fins: 3 per stage, as in the model). Variants: 4 fins per stage, and
fin planforms scaled by 0.9 and 1.1 (planform uncertainty of the model
measurement). For each variant the booster is recalibrated on V01 (Isp, t_b1),
then the held-out VSB-30 apogee errors, the unit-wind effect and the dwell and
burnout Mach number of the depressed Concepts A and D are computed, as in
run_aero_sensitivity.py.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import calibration, flight, geometry as G, stability as St, vehicle  # noqa: E402
import run_aero_sensitivity as ra  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
MACHS = [0.3, 0.6, 0.9, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]


def tables(n_fins=3, fin_scale=1.0):
    g = G.vsb30_geometry(n_fins, fin_scale=fin_scale)
    veh = vehicle.vsb30(400.0)
    items = G.vsb30_mass_items(veh, g["stations"])
    s31, s30 = veh.stages
    out = {}
    for name, cfg, t in (("full, lift-off", g["full"], 0.0), ("full, S31 burnout", g["full"], s31.t_drop - 1e-3),
                         ("upper, S30 ignition", g["upper"], s30.t_ign + 1e-3),
                         ("upper, S30 burnout", g["upper"], s30.t_ign + s30.burn_time)):
        m, xcg, I = G.mass_properties(items, t)
        rows = []
        for M in MACHS:
            d = St.pitch_derivatives(cfg, M, xcg)
            rows.append({"mach": M, "cn_alpha": d["cn_alpha"], "x_cp": d["x_cp"], "margin_cal": d["margin_cal"],
                         "cm_alpha": d["cm_alpha"], "cmq": d["cmq"]})
        out[name] = {"t": t, "mass": m, "x_cg": xcg, "inertia": I, "rows": rows}
    geo = {"stations": g["stations"], "fins": [{k: getattr(f, k) for k in ("name", "x_le", "c_root", "c_tip", "span",
                                                                          "x_tip_offset", "n")}
                                              for f in g["full"].fins]}
    return out, geo


def mission(factory):
    flight.AERO_FACTORY = factory
    try:
        c = calibration.calibrate(calibration.V01)
        r = {"isp": c["isp"], "b1": c["b1"], "uwe_km_per_mps": ra.uwe(None, c["isp"], c["b1"]),
             "held_out": ra.held_out(c["isp"], c["b1"])}
        r.update(ra.concepts(c["isp"], c["b1"]))
        return r
    finally:
        flight.AERO_FACTORY = None


if __name__ == "__main__":
    t0 = time.time()
    tab, geo = tables()
    out = {"geometry": geo, "tables": tab, "variants": {}}
    out["tables_4fins"] = tables(4)[0]
    variants = {"generic": None,
                "geometric, 3 fins": lambda v: St.vsb30_aero(v, 3),
                "geometric, 4 fins": lambda v: St.vsb30_aero(v, 4),
                "geometric, fins x0.9": lambda v: St.vsb30_aero(v, 3, fin_scale=0.9),
                "geometric, fins x1.1": lambda v: St.vsb30_aero(v, 3, fin_scale=1.1)}
    for name, fac in variants.items():
        out["variants"][name] = mission(fac)
        print(name, json.dumps(out["variants"][name]), flush=True)
    geo_keys = [k for k in variants if k != "generic"]
    for cc in "AD":
        d = [out["variants"][k]["dwell_" + cc] for k in geo_keys]
        out["dwell_range_" + cc] = [min(d), max(d)]
    out["cpu_s"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "geometric_aero.json"), "w") as f:
        json.dump(out, f, indent=1)
