"""Validation domain: flight envelope of the ten VSB-30 validation flights
(launch elevation, payload, maximum ascent dynamic pressure, maximum q*alpha
during the boost, burnout Mach number, apogee) compared with the nominal
concepts of the case study (results/casestudy.json).

Elevation, payload and apogee of the flights are flight data (the measured
apogee of Garcia et al. 2011); q, q-alpha and burnout Mach number are SHyMA
reconstructions of the flight configurations. The reconstructed apogee is
kept as apogee_sim_km for reference."""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import earth, flight, vehicle  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def envelope(tr, t_burnout):
    ta = tr.apogee()[0]
    asc = [k for k in range(len(tr.t)) if tr.t[k] <= ta]
    na = min(len(tr.alpha), len(tr.t)) if getattr(tr, "alpha", None) else 0
    qa = max((tr.q[k] * tr.alpha[k] for k in range(na) if tr.thrust[k] > 0 or tr.t[k] < 60.0), default=0.0)
    kb = min(range(len(tr.t)), key=lambda k: abs(tr.t[k] - t_burnout))
    return {"q_max_ascent_kpa": max(tr.q[k] for k in asc) / 1e3, "qalpha_max_kpa_deg": qa / 1e3,
            "burnout_mach": tr.mach[kb], "apogee_km": tr.apogee()[1] / 1e3}


if __name__ == "__main__":
    rows = [r for r in csv.DictReader(l for l in open(os.path.join(ROOT, "data", "vsb30_flights_garcia2011.csv"))
                                      if not l.startswith("#"))]
    fl = []
    for r in rows:
        v = vehicle.vsb30(float(r["payload_kg"]))
        tr = flight.fly(v, earth.STATIONS[r["site"]], float(r["elevation_deg"]), float(r["azimuth_deg"]), 60.0,
                        record_every=1)
        e = envelope(tr, v.t_burnout)
        e["apogee_sim_km"] = e["apogee_km"]
        e["apogee_km"] = float(r["apogee_actual_km"])
        e.update({"flight": r["flight"], "site": r["site"], "elevation": float(r["elevation_deg"]),
                  "payload_kg": float(r["payload_kg"])})
        fl.append(e)
        print(e, flush=True)
    cs = json.load(open(os.path.join(ROOT, "results", "casestudy.json")))
    conc = {c: {"elevation": n["elevation"], "payload_kg": n["payload"], "q_max_ascent_kpa": n["q_max_ascent_kpa"],
                "qalpha_max_kpa_deg": n.get("qalpha_max_kpa_deg"), "burnout_mach": n["burnout_mach"],
                "apogee_km": n["apogee_km"]} for c, n in cs["nominal"].items()}
    keys = ("elevation", "payload_kg", "q_max_ascent_kpa", "qalpha_max_kpa_deg", "burnout_mach", "apogee_km")
    rng = {k: [min(f[k] for f in fl), max(f[k] for f in fl)] for k in keys}
    inside = {c: {k: rng[k][0] - 1e-9 <= conc[c][k] <= rng[k][1] + 1e-9 for k in keys} for c in conc}
    with open(os.path.join(ROOT, "results", "validation_domain.json"), "w") as f:
        json.dump({"flights": fl, "flight_range": rng, "concepts": conc, "inside": inside}, f, indent=1)
    print(json.dumps(rng), json.dumps(inside))
