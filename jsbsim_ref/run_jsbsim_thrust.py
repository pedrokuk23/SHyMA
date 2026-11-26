"""Code-to-code verification with propulsion: SHyMA versus JSBSim.

A point-mass body (400 kg empty, 600 kg propellant, C_D S = 0.05 m^2) is launched
vertically from 1 km above the Alcantara Launch Center with 100 m/s and burns for
30 s with a constant mass flow (20 kg/s), a vacuum specific impulse of 250 s and an
ambient-pressure thrust loss over a 0.10 m^2 exit area, F = mdot Isp g0 - p_a A_e,
then coasts. In JSBSim the thrust is an external body-axis force and the propellant a
point mass, both updated every step (1 ms) from the JSBSim atmosphere; the body
attitude is held along the launch vertical (no aerodynamic moments). SHyMA flies the
same case in its point-mass mode (thrust along the air-relative velocity). The two
trajectories are compared at equal altitude (Mach number, dynamic pressure, time)
during the powered ascent and the coast, up to 120 km.
"""
import json
import math
import os
import sys

import jsbsim

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from shyma import earth, trajectory, vehicle  # noqa: E402
from shyma.atmosphere import atmosphere  # noqa: E402

LAT, LON = earth.STATIONS["CLA"][:2]
FT, LBF, G0 = 0.3048, 4.4482216152605, 9.80665
M_DRY, M_PROP, MDOT, ISP, AE, CDS = 400.0, 600.0, 20.0, 250.0, 0.10, 0.05
TB = M_PROP / MDOT
H0, V0 = 1000.0, 100.0


def run_jsbsim(dt=0.001, h_end=120e3):
    fdm = jsbsim.FGFDMExec(None)
    fdm.set_debug_level(0)
    fdm.set_aircraft_path(os.path.join(HERE, "aircraft"))
    fdm.load_model("thruster")
    fdm.set_dt(dt)
    fdm["ic/lat-geod-deg"] = LAT
    fdm["ic/long-gc-deg"] = LON
    fdm["ic/h-sl-ft"] = H0 / FT
    fdm["ic/vn-fps"] = 0.0
    fdm["ic/ve-fps"] = 0.0
    fdm["ic/vd-fps"] = -V0 / FT
    fdm["ic/theta-deg"] = 90.0
    fdm["ic/psi-true-deg"] = 0.0
    fdm["ic/phi-deg"] = 0.0
    fdm.run_ic()
    out = {"t": [], "alt_km": [], "mach": [], "q_kpa": []}
    k = 0
    while fdm["position/h-sl-ft"] * FT < h_end and fdm.get_sim_time() < 400:
        t = fdm.get_sim_time()
        if t < TB:
            p = fdm["atmosphere/P-psf"] * 47.880259
            fdm["external_reactions/thrust/magnitude"] = max(MDOT * ISP * G0 - p * AE, 0.0) / LBF
            fdm["inertia/pointmass-weight-lbs[0]"] = (M_PROP - MDOT * (t + 0.5 * dt)) / 0.45359237
        else:
            fdm["external_reactions/thrust/magnitude"] = 0.0
            fdm["inertia/pointmass-weight-lbs[0]"] = 0.0
        if k % 20 == 0:
            out["t"].append(t)
            out["alt_km"].append(fdm["position/h-sl-ft"] * FT / 1e3)
            out["mach"].append(fdm["velocities/mach"])
            out["q_kpa"].append(fdm["aero/qbar-psf"] * 47.880259 / 1e3)
        fdm.run()
        k += 1
        if fdm["velocities/v-down-fps"] > 0 and t > TB:      # apogee passed
            break
    return out


class Body(vehicle.Vehicle):
    def drag_area(self, t, mach, tc=None, re_per_m=None):
        return CDS


def run_shyma(h_end=120e3):
    st = vehicle.Stage("S", M_PROP, 0.0, TB, ISP, t_ign=0.0, t_drop=1e9, exit_area=AE)
    veh = Body([st], M_DRY, 1.0)
    x, y, z = earth.llh_to_ecef(LAT, LON, H0)
    e, n, u = earth.enu_basis(LAT, LON)
    s0 = [x, y, z, V0 * u[0] - earth.OMEGA * y, V0 * u[1] + earth.OMEGA * x, V0 * u[2]]
    tr = trajectory.propagate(veh, 0.0, s0, dt_atm=0.001, dt_vac=0.01, record_every=20,
                              t_end=160.0)
    return tr


def at_alt(ts, hs, ys, h):
    for k in range(len(hs) - 1):
        if hs[k] <= h <= hs[k + 1] and hs[k + 1] > hs[k]:
            w = (h - hs[k]) / (hs[k + 1] - hs[k])
            return ys[k] + w * (ys[k + 1] - ys[k])
    return None


if __name__ == "__main__":
    js = run_jsbsim()
    tr = run_shyma()
    sh = {"t": tr.t, "alt_km": [a / 1e3 for a in tr.alt], "mach": tr.mach, "q_kpa": [q / 1e3 for q in tr.q]}
    rows = []
    for h in [2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]:
        r = {"alt_km": h}
        for key in ("mach", "q_kpa", "t"):
            a, b = at_alt(sh["t"], sh["alt_km"], sh[key], h), at_alt(js["t"], js["alt_km"], js[key], h)
            r[key] = (a, b)
        rows.append(r)
    dm = max(abs(r["mach"][0] / r["mach"][1] - 1) for r in rows if r["mach"][1])
    dq = max(abs(r["q_kpa"][0] / r["q_kpa"][1] - 1) for r in rows if r["q_kpa"][1] and r["q_kpa"][1] > 1e-3)
    dt_ = max(abs(r["t"][0] - r["t"][1]) for r in rows)
    apo_s = max(sh["alt_km"])
    out = {"case": {"m_dry": M_DRY, "m_prop": M_PROP, "mdot": MDOT, "isp_vac": ISP, "exit_area_m2": AE,
                    "cds_m2": CDS, "h0_m": H0, "v0_mps": V0, "burn_s": TB},
           "max_rel_mach": dm, "max_rel_q": dq, "max_dt_s": dt_, "rows": rows,
           "burnout_alt_km": at_alt(sh["alt_km"], sh["t"], sh["alt_km"], TB) if False else None,
           "max_alt_compared_km": min(apo_s, max(js["alt_km"]))}
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    for r in rows:
        print(r)
    json.dump(out, open(os.path.join(HERE, "..", "results", "jsbsim_thrust.json"), "w"), indent=1)
