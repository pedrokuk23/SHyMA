"""Sensitivity of the mission results to the generic aerodynamic representation
of the booster: stability coefficients (C_N_alpha, static margin, C_mq),
each halved and doubled, and the zero-lift drag curve of the attached stack,
scaled by 0.8 and 1.2.

Each variant is evaluated twice:
  * frozen: propulsion kept at the nominal V01 calibration (direct effect of
    the aerodynamic change);
  * recalibrated: Isp and t_b1 re-fitted on V01 for the variant (the effect
    that remains after the calibration the method allows, which can absorb
    part of the aerodynamic change).
For both, the unit-wind effect (impact displacement for a 1 m/s wind from
60 deg in the lowest 3 km, V01), the apogee error on the nine held-out
VSB-30 flights and the dwell and burnout Mach number of the depressed
Concepts A and D (largest q-alpha, outside the elevation range of the
validation flights) are computed.
"""
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from shyma import calibration, earth, flight, vehicle, wind  # noqa: E402
import run_casestudy as cs  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
BASE = dict(flight.VSB30_AERO)
SITE = earth.STATIONS["CLA"]
CD_SCALE = [1.0]
_drag_area = vehicle.Vehicle.drag_area


def _scaled_drag_area(self, t, mach, tc=None, re_per_m=None):
    a = _drag_area(self, t, mach, tc, re_per_m)
    return a * CD_SCALE[0] if self.attached(t if tc is None else tc) else a


vehicle.Vehicle.drag_area = _scaled_drag_area


def held_out(isp, b1):
    rows = [r for r in csv.DictReader(l for l in open(os.path.join(ROOT, "data", "vsb30_flights_garcia2011.csv"))
                                      if not l.startswith("#")) if r["flight"] != "V01"]
    err = []
    for r in rows:
        v = vehicle.vsb30(float(r["payload_kg"]), isp=isp, s31_burn=b1)
        tr = flight.fly(v, earth.STATIONS[r["site"]], float(r["elevation_deg"]), float(r["azimuth_deg"]), 60.0,
                        record_every=50)
        err.append(100 * (tr.apogee()[1] / 1e3 / float(r["apogee_actual_km"]) - 1))
    return {"mean_pct": sum(err) / len(err), "rms_pct": math.sqrt(sum(e * e for e in err) / len(err))}


def concepts(isp, b1):
    out = {}
    for cc in "AD":
        v, tr = cs.fly(cc, 400.0, cs.NOMINAL[cc], isp=isp, s31_burn=b1, rec=1)
        sm = cs.summarize(cc, v, tr, light=True)
        out["dwell_" + cc] = sm["window_s"]
        out["burnout_mach_" + cc] = sm["burnout_mach"]
    return out


def uwe(aero, isp, b1):
    w1 = wind.uniform_wind(1.0, 60.0, top_m=3000.0)
    out = []
    for w in (None, w1):
        v = vehicle.vsb30(392.0, isp=isp, s31_burn=b1)
        tr = flight.fly(v, SITE, 84.0, 45.0, 60.0, wind=w, aero=aero, record_every=50)
        out.append((tr.lat[-1], tr.lon[-1]))
    return earth.great_circle_km(out[0][0], out[0][1], out[1][0], out[1][1])


if __name__ == "__main__":
    res = {}
    variants = [(None, 1.0), ("cn_alpha", 0.5), ("cn_alpha", 2.0), ("margin_cal", 0.5), ("margin_cal", 2.0),
                ("cmq", 0.5), ("cmq", 2.0), ("cd", 0.8), ("cd", 1.2)]
    nom = None
    for key, f in variants:
        aero = dict(BASE)
        if key in BASE:
            aero[key] = BASE[key] * f
        CD_SCALE[0] = f if key == "cd" else 1.0
        flight.VSB30_AERO.clear()
        flight.VSB30_AERO.update(aero)
        name = "nominal" if key is None else "%s_x%.1f" % (key, f)
        c = calibration.calibrate(calibration.V01)
        r = {"isp": c["isp"], "b1": c["b1"], "uwe_km_per_mps": uwe(aero, c["isp"], c["b1"]),
             "held_out": held_out(c["isp"], c["b1"])}
        r.update(concepts(c["isp"], c["b1"]))
        if key == "cd":      # launch settings re-tuned for the variant (grids of run_thrust_shape.py)
            ra = max(((el, cs.summarize("A", *cs.fly("A", 400.0, el, isp=c["isp"], s31_burn=c["b1"], rec=1),
                                        light=True)["window_s"]) for el in [56 + 0.5 * i for i in range(25)]),
                     key=lambda z: z[1])
            rd = max(((el, co, cs.summarize("D", *cs.fly("D", 400.0, el, coast=co, isp=c["isp"], s31_burn=c["b1"],
                                                         rec=1), light=True)["window_s"])
                      for el in (69.0, 71.0, 73.0, 75.0, 77.0) for co in [17 + i for i in range(14)]),
                     key=lambda z: z[2])
            r["A_retuned"] = {"el": ra[0], "dwell_s": ra[1]}
            r["D_retuned"] = {"el": rd[0], "coast": rd[1], "dwell_s": rd[2]}
        if nom is None:
            nom = (c["isp"], c["b1"])
        else:
            fr = {"held_out": held_out(*nom)}
            fr.update(concepts(*nom))
            r["frozen"] = fr
        res[name] = r
        print(name, json.dumps(r), flush=True)
    CD_SCALE[0] = 1.0
    flight.VSB30_AERO.clear()
    flight.VSB30_AERO.update(BASE)
    n = res["nominal"]["uwe_km_per_mps"]
    stab = [k for k in res if k != "nominal" and not k.startswith("cd_")]
    ch = [100 * (res[k]["uwe_km_per_mps"] / n - 1) for k in stab]
    res["uwe_change_pct_range"] = [min(ch), max(ch)]
    for cc in "AD":
        d = [res[k]["dwell_" + cc] for k in stab]
        res["dwell_range_" + cc] = [min(d), max(d)]
        for grp, keys in (("stab", stab), ("cd", ["cd_x0.8", "cd_x1.2"])):
            for mode in ("recal", "frozen"):
                d = [(res[k] if mode == "recal" else res[k]["frozen"])["dwell_" + cc] for k in keys]
                res["dwell_%s_%s_%s" % (grp, mode, cc)] = [min(d), max(d)]
    for grp, keys in (("stab", stab), ("cd", ["cd_x0.8", "cd_x1.2"])):
        for mode in ("recal", "frozen"):
            d = [abs((res[k] if mode == "recal" else res[k]["frozen"])["held_out"]["mean_pct"]) for k in keys]
            res["held_out_absmean_max_%s_%s" % (grp, mode)] = max(d)
    for cc in "AD":
        d = [res[k][cc + "_retuned"]["dwell_s"] for k in ("cd_x0.8", "cd_x1.2")]
        res["dwell_cd_retuned_" + cc] = [min(d), max(d)]
    with open(os.path.join(ROOT, "results", "aero_sensitivity.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(res["uwe_change_pct_range"])
