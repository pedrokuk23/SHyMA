"""Code-to-code verification of the SHyMA trajectory core against JSBSim
(open-source flight dynamics engine that took part in the NASA NESC 6-DOF
check-case study, NASA/TM-2015-218675). Case: ballistic body, m = 500 kg,
Cd*S = 0.1 m^2 (beta = 5000 kg/m^2), entering at 85 km above Alcantara
with 2.4 km/s Earth-relative speed, -70 deg flight-path angle, heading 45 deg.
"""
import json
import math
import os
import sys

import jsbsim

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from shyma import analysis, earth, trajectory, vehicle  # noqa: E402

LAT, LON = earth.STATIONS["CLA"][:2]
CASES = {"entry85": (85e3, 2400.0, -70.0, 45.0), "apogee250": (250e3, 700.0, 0.0, 45.0)}
H0, V0, G0, PSI = CASES["entry85"]
FT = 0.3048


def run_jsbsim(dt=0.002):
    fdm = jsbsim.FGFDMExec(None)
    fdm.set_debug_level(0)
    fdm.set_aircraft_path(os.path.join(HERE, "aircraft"))
    fdm.set_engine_path(os.path.join(HERE, "aircraft"))
    fdm.set_systems_path(os.path.join(HERE, "aircraft"))
    fdm.load_model("sphere")
    fdm.set_dt(dt)
    vn = V0 * math.cos(math.radians(G0)) * math.cos(math.radians(PSI))
    ve = V0 * math.cos(math.radians(G0)) * math.sin(math.radians(PSI))
    vd = -V0 * math.sin(math.radians(G0))
    fdm["ic/lat-geod-deg"] = LAT
    fdm["ic/long-gc-deg"] = LON
    fdm["ic/h-sl-ft"] = H0 / FT
    fdm["ic/vn-fps"] = vn / FT
    fdm["ic/ve-fps"] = ve / FT
    fdm["ic/vd-fps"] = vd / FT
    fdm["ic/psi-true-deg"] = PSI
    fdm["ic/theta-deg"] = G0
    fdm["ic/phi-deg"] = 0.0
    fdm.run_ic()
    out = {"t": [], "alt_km": [], "mach": [], "q_kpa": [], "lat": [], "lon": []}
    k = 0
    while fdm["position/h-sl-ft"] * FT > 0.0 and fdm.get_sim_time() < 1000:
        if k % 50 == 0:
            out["t"].append(fdm.get_sim_time())
            out["alt_km"].append(fdm["position/h-sl-ft"] * FT / 1e3)
            out["mach"].append(fdm["velocities/mach"])
            out["q_kpa"].append(fdm["aero/qbar-psf"] * 47.880259 / 1e3)
            out["lat"].append(fdm["position/lat-geod-deg"])
            out["lon"].append(fdm["position/long-gc-deg"])
        fdm.run()
        k += 1
    out["t"].append(fdm.get_sim_time())
    out["alt_km"].append(fdm["position/h-sl-ft"] * FT / 1e3)
    out["mach"].append(fdm["velocities/mach"])
    out["q_kpa"].append(fdm["aero/qbar-psf"] * 47.880259 / 1e3)
    out["lat"].append(fdm["position/lat-geod-deg"])
    out["lon"].append(fdm["position/long-gc-deg"])
    return out


def run_shyma():
    body = vehicle.BallisticBody(500.0, 5000.0)
    # Earth-relative initial velocity -> inertial state
    x, y, z = earth.llh_to_ecef(LAT, LON, H0)
    e, n, u = earth.enu_basis(LAT, LON)
    g, p = math.radians(G0), math.radians(PSI)
    d = [math.cos(g) * math.sin(p) * e[i] + math.cos(g) * math.cos(p) * n[i] + math.sin(g) * u[i]
         for i in range(3)]
    s0 = [x, y, z, V0 * d[0] - earth.OMEGA * y, V0 * d[1] + earth.OMEGA * x, V0 * d[2]]
    tr = trajectory.propagate(body, 0.0, s0, dt_atm=0.01)
    return tr


def interp(xs, ys, x):
    for k in range(len(xs) - 1):
        if (xs[k] - x) * (xs[k + 1] - x) <= 0 and xs[k] != xs[k + 1]:
            w = (x - xs[k]) / (xs[k + 1] - xs[k])
            return ys[k] + w * (ys[k + 1] - ys[k])
    return float("nan")


def compare(name):
    global H0, V0, G0, PSI
    H0, V0, G0, PSI = CASES[name]
    js = run_jsbsim()
    tr = run_shyma()
    am = {"t": tr.t, "alt_km": [a / 1e3 for a in tr.alt], "mach": tr.mach,
          "q_kpa": [q / 1e3 for q in tr.q]}
    cmp = []
    for h in (80, 70, 60, 50, 40, 35, 30, 25, 20, 15, 10, 5):
        r = {"alt_km": h}
        for src, d in (("jsbsim", js), ("shyma", am)):
            r[src] = {f: interp(d["alt_km"], d[f], h) for f in ("t", "mach", "q_kpa")}
        r["d_mach_pct"] = 100 * (r["shyma"]["mach"] - r["jsbsim"]["mach"]) / r["jsbsim"]["mach"]
        r["d_q_pct"] = 100 * (r["shyma"]["q_kpa"] - r["jsbsim"]["q_kpa"]) / r["jsbsim"]["q_kpa"]
        r["d_t_s"] = r["shyma"]["t"] - r["jsbsim"]["t"]
        cmp.append(r)
        print(h, "M %.3f %.3f (%.2f%%)  q %.2f %.2f (%.2f%%)  t %.2f %.2f" % (
            r["jsbsim"]["mach"], r["shyma"]["mach"], r["d_mach_pct"],
            r["jsbsim"]["q_kpa"], r["shyma"]["q_kpa"], r["d_q_pct"], r["jsbsim"]["t"], r["shyma"]["t"]))
    imp = earth.great_circle_km(js["lat"][-1], js["lon"][-1], tr.lat[-1], tr.lon[-1])
    print(name, "impact separation km", imp, "t_impact", js["t"][-1], tr.t[-1])
    return {"case": {"h0_km": H0 / 1e3, "v0": V0, "fpa0": G0, "heading": PSI, "beta": 5000.0},
           "comparison": cmp, "impact_sep_km": imp, "t_impact_jsbsim": js["t"][-1],
           "t_impact_shyma": tr.t[-1], "jsbsim_version": jsbsim.__version__,
           "jsbsim_profile": js, "shyma_profile": {k: v[::20] for k, v in am.items()}}


if __name__ == "__main__":
    earth.set_model("sphere")   # spherical geodesy of the earlier ASMAT satellite tool
    res = {n: compare(n) for n in CASES}
    # same cases with the WGS-84 geodesy option (JSBSim's own Earth model)
    earth.set_model("wgs84")
    for n in list(CASES):
        res[n + "_wgs84"] = compare(n)
    earth.reset_model()
    with open(os.path.join(HERE, "..", "results", "jsbsim_verification.json"), "w") as f:
        json.dump(res, f)
