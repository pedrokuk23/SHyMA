"""Verification against the NASA NESC atmospheric check-cases
(Murri, Jackson & Shelton, NASA/TM-2015-218675; data from the NESC Academy).

Case 05: sphere dropped from 30,000 ft, round rotating Earth, inverse-square gravity.
Case 09: sphere launched eastward from (0 N, 0 E), WGS-84, J2 gravity.
Case 10: sphere launched northward from (0 N, 0 E), WGS-84, J2 gravity.
Case 17: unguided two-stage rocket, WGS-84, J2; compared from t = 20 s
         (after the 6-DOF lift-off weathercocking transient) up to second-stage
         ignition at t = 134.17 s, i.e., over the atmospheric boost, burnout and
         the 97 s exo-atmospheric coast. Thrust is along the air-relative
         velocity (zero angle of attack); after staging the 6-DOF vehicle keeps
         an inertially fixed attitude, which a point-mass model cannot represent.

The sphere has m = 1 slug, S = 0.1963495 ft^2, CD = 0.1 (Vol. II, Table 2).
The two-stage rocket data (Vol. II, Sec. B.1.6): stage 1 thrust 17 MN,
Isp 360 s, 180 t of propellant, 35 t dry; stage 2 thrust 5 MN, Isp 390 s,
80 t of propellant, 15 t dry; 4 t payload; 96.79 s coast; reference area 7 m^2;
US76 atmosphere with zero density above 86 km. The zero-angle-of-attack drag
coefficient (0.212) is taken from the reference aerodynamic force output,
because the drag table itself is distributed only in the DAVE-ML model files.
"""
import csv
import glob
import json
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import earth, trajectory, vehicle  # noqa: E402
from shyma.atmosphere import atmosphere  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA = os.path.join(ROOT, "validation_data")
FT = 0.3048
SLUG = 14.59390294

SPHERE_M = 1.0 * SLUG
SPHERE_CDA = 0.1 * 0.1963495 * FT ** 2
SPHERE_BETA = SPHERE_M / SPHERE_CDA


def load_case(pattern):
    sims = {}
    for f in sorted(glob.glob(os.path.join(DATA, pattern))):
        name = os.path.basename(f).split("_sim_")[1].split(".")[0]
        with open(f) as fh:
            rows = list(csv.DictReader(fh))
        sims["SIM " + name] = rows
    return sims


def col(rows, key, scale=1.0):
    return [float(r[key]) * scale for r in rows]


def interp(ts, ys, t):
    if t <= ts[0]:
        return ys[0]
    if t >= ts[-1]:
        return ys[-1]
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t:
            lo = mid
        else:
            hi = mid
    w = (t - ts[lo]) / (ts[hi] - ts[lo])
    return ys[lo] + w * (ys[hi] - ys[lo])


VARS = {  # name: (csv column, scale to SI, trajectory attribute, scale)
    "alt_m": ("altitudeMsl_ft", FT, "alt", 1.0),
    "vn_mps": ("feVelocity_ft_s_X", FT, "vn", 1.0),
    "ve_mps": ("feVelocity_ft_s_Y", FT, "ve", 1.0),
    "vd_mps": ("feVelocity_ft_s_Z", FT, "vd", 1.0),
    "lat_deg": ("latitude_deg", 1.0, "lat", 1.0),
    "lon_deg": ("longitude_deg", 1.0, "lon", 1.0),
    "mach": ("mach", 1.0, "mach", 1.0),
}


def series(sims, c, sc):
    """Time series of column c for every simulation that provides it."""
    out = []
    for rows in sims.values():
        if c in rows[0] and rows[0][c] not in ("", None):
            out.append((col(rows, "time"), col(rows, c, sc)))
    return out


def compare(sims, tr, times, tmax=None):
    out = {}
    for v, (c, sc, attr, sa) in VARS.items():
        ser = series(sims, c, sc)
        if len(ser) < 2:
            continue
        errs, spreads = [], []
        for t in times:
            if tmax and v in tmax and t > tmax[v]:
                continue
            refs = [interp(ts, ys, t) for ts, ys in ser]
            med = statistics.median(refs)
            mine = interp(tr.t, [x * sa for x in getattr(tr, attr)], t)
            errs.append(mine - med)
            spreads.append(max(refs) - min(refs))
        out[v] = {"n_ref": len(ser), "max_abs_err": max(abs(e) for e in errs), "max_spread": max(spreads),
                  "final_err": errs[-1], "final_spread": spreads[-1]}
    return out


def ground_distance_m(sims, tr, t):
    """Distance between SHyMA and the median reference position at time t."""
    lats = [interp(col(r, "time"), col(r, "latitude_deg"), t) for r in sims.values()]
    lons = [interp(col(r, "time"), col(r, "longitude_deg"), t) for r in sims.values()]
    la, lo = statistics.median(lats), statistics.median(lons)
    return earth.great_circle_km(la, lo, interp(tr.t, tr.lat, t), interp(tr.t, tr.lon, t)) * 1e3


def case05():
    earth.set_model("round", j2=False, rotation=True)
    body = vehicle.BallisticBody(SPHERE_M, SPHERE_BETA)
    s0 = trajectory.state_from_ned(0.0, 0.0, 30000 * FT, 0.0, 0.0, 0.0)
    tr = trajectory.propagate(body, 0.0, s0, dt_atm=0.01, t_end=30.0, stop_alt=-1e9)
    earth.reset_model()
    return tr


def case_cannon(vn, ve):
    earth.set_model("wgs84")
    body = vehicle.BallisticBody(SPHERE_M, SPHERE_BETA)
    s0 = trajectory.state_from_ned(0.0, 0.0, 0.0, vn * FT, ve * FT, -1000 * FT)
    tr = trajectory.propagate(body, 0.0, s0, dt_atm=0.01, t_end=30.0, stop_alt=-1e9)
    earth.reset_model()
    return tr


class TwoStage(vehicle.Vehicle):
    CD0 = 0.212

    def drag_area(self, t, mach, tc=None):
        return self.CD0 * self.ref_area


def atm_us76_cut86(alt_km, level="mean", rho_factor=1.0):
    T, p, rho, a = atmosphere(min(alt_km, 86.0), level, rho_factor)
    if alt_km > 86.0:
        return T, 0.0, 0.0, a
    return T, p, rho, a


def case17(t0=20.0, t_end=134.17, init="median"):
    """Two-stage rocket from t0. init = 'median': component-wise median of the reference ECI
    states at t0 (all simulations); otherwise the name of one simulation (e.g. 'SIM 05')."""
    earth.set_model("wgs84")
    g0 = 9.80665
    md1, md2 = 17e6 / (360 * g0), 5e6 / (390 * g0)
    b1, b2 = 180000.0 / md1, 80000.0 / md2
    t_stage = b1 + 96.79
    s1 = vehicle.Stage("S1", 180000.0, 35000.0, b1, 360.0, 0.0, t_stage)
    s2 = vehicle.Stage("S2", 80000.0, 15000.0, b2, 390.0, t_stage, 1e9)
    veh = TwoStage([s1, s2], 4000.0, 7.0)
    sims = load_case("Atmos_17*/*.csv")
    states = []
    for name, ref in sims.items():
        if init != "median" and name != init:
            continue
        r0 = min(ref, key=lambda r: abs(float(r["time"]) - t0))
        states.append([float(r0["eiPosition_ft_" + k]) * FT for k in "XYZ"] +
                      [float(r0["eiVelocity_ft_s_" + k]) * FT for k in "XYZ"])
    s0 = [statistics.median(c) for c in zip(*states)]
    tr = trajectory.propagate(veh, t0, s0, dt_atm=0.001, dt_vac=0.01,
                              t_end=t_end, stop_alt=-1e9, atm=atm_us76_cut86, record_every=10)
    earth.reset_model()
    return tr, {"stage1_burn_s": b1, "t_staging_s": t_stage, "stage2_burn_s": b2}


if __name__ == "__main__":
    res = {}
    for name, fn, pat in (("case05", case05, "Atmos_05*/*.csv"),
                          ("case09", lambda: case_cannon(0.0, 1000.0), "Atmos_09*/*.csv"),
                          ("case10", lambda: case_cannon(1000.0, 0.0), "Atmos_10*/*.csv")):
        sims = load_case(pat)
        tr = fn()
        times = [0.5 * k for k in range(1, 61)]
        res[name] = {"n_sims": len(sims), "sims": list(sims), "metrics": compare(sims, tr, times),
                     "profile": {"t": tr.t[::10], "alt_m": tr.alt[::10]},
                     "ref_profiles": {k: {"t": col(r, "time")[::5], "alt_m": col(r, "altitudeMsl_ft", FT)[::5]}
                                      for k, r in sims.items()},
                     "pos_err_final_m": ground_distance_m(sims, tr, 30.0) if name != "case05" else None}
        print(name, json.dumps(res[name]["metrics"], indent=0)[:900], res[name]["pos_err_final_m"])
    sims = load_case("Atmos_17*/*.csv")
    tr, info = case17()
    times = [20.0 + 2.0 * k for k in range(0, 58)] + [134.17]
    res["case17"] = {"n_sims": len(sims), "sims": list(sims), "info": info,
                     "metrics": compare(sims, tr, times, tmax={"mach": 75.0}),
                     "pos_err_staging_m": ground_distance_m(sims, tr, 134.17),
                     "profile": {"t": tr.t, "alt_m": tr.alt, "v": [math.sqrt(a * a + b * b + c * c)
                                                                 for a, b, c in zip(tr.vn, tr.ve, tr.vd)]},
                     "ref_profiles": {k: {"t": col(r, "time")[::100], "alt_m": col(r, "altitudeMsl_ft", FT)[::100],
                                          "v": [math.sqrt(a * a + b * b + c * c) for a, b, c in zip(
                                              col(r, "feVelocity_ft_s_X", FT)[::100], col(r, "feVelocity_ft_s_Y", FT)[::100],
                                              col(r, "feVelocity_ft_s_Z", FT)[::100])]}
                                      for k, r in sims.items()},
                     "checkpoints": {}}
    # sensitivity to the initial state: start from each simulation in turn
    res["case17"]["init_sensitivity"] = {}
    for name in sims:
        tri, _ = case17(init=name)
        res["case17"]["init_sensitivity"][name] = compare(sims, tri, times, tmax={"mach": 75.0})
    for tc in (37.38, 80.0, 134.17):
        cp = {}
        for v, (c, sc, attr, sa) in VARS.items():
            refs = [interp(ts, ys, tc) for ts, ys in series(sims, c, sc)]
            if not refs:
                continue
            cp[v] = {"ref_median": statistics.median(refs), "ref_min": min(refs), "ref_max": max(refs),
                     "shyma": interp(tr.t, getattr(tr, attr), tc)}
        res["case17"]["checkpoints"]["%.2f" % tc] = cp
    print("case17", json.dumps(res["case17"]["metrics"], indent=0)[:1200], res["case17"]["pos_err_staging_m"])
    print(json.dumps(res["case17"]["checkpoints"], indent=0)[:2500])
    with open(os.path.join(ROOT, "results", "nesc_checkcases.json"), "w") as f:
        json.dump(res, f, indent=1)
