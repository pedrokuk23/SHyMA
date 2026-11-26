"""Code verification of the pseudo-5-DOF attitude dynamics against the linear
short-period (weathercock) solution.

With gravity, Earth rotation, thrust and axial drag switched off and a
constant atmosphere, a small initial angle of attack alpha0 of a coasting
rocket obeys

    alpha'' + (d + n) alpha' + (k + d n) alpha = 0,
    k = l q S CNa / I,   d = -q S d_ref^2 Cmq / (2 V I),   n = q S CNa / (m V),

with alpha'(0) = -n alpha0 (the body starts at rest). The integrator of
shyma.trajectory5 (full 12-state RK4) is compared with this solution for the
generic coefficients of the paper and for a full stack at Mach 3 at 20 km.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import earth, trajectory5 as t5, vehicle  # noqa: E402
from shyma.trajectory import Env  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def run(alpha0_deg=1.0, h_km=20.0, mach=3.0, dt=1e-3, t_end=8.0):
    T, p, rho, a = 216.65, 5474.9, 0.088035, 295.07          # US76 at 20 km
    env = Env(atm=lambda alt_km, level, rf: (T, p, rho, a))
    veh = vehicle.vsb30(400.0, cd_scale=0.0)
    for s in veh.stages:                                     # coasting, full stack, no thrust
        s.t_ign, s.burn_time, s.t_drop = 1e9, 1.0, 2e9
    aero = t5.Aero5(0.557)
    g_saved, om_saved = earth.gravity, earth.OMEGA
    earth.gravity = lambda x, y, z: (0.0, 0.0, 0.0)
    earth.OMEGA = 0.0
    try:
        r0 = earth.RE + h_km * 1e3
        V = mach * a
        a0 = math.radians(alpha0_deg)
        s = [r0, 0.0, 0.0, 0.0, V, 0.0, 0.0, math.cos(a0), math.sin(a0), 0.0, 0.0, 0.0]
        m = veh.mass(0.0)
        I = aero.inertia(veh, 0.0)
        S = veh.ref_area
        q = 0.5 * rho * V * V
        k = aero.margin * q * S * aero.cn_alpha / I
        d = -q * S * aero.d ** 2 * aero.cmq / (2 * V * I)
        n = q * S * aero.cn_alpha / (m * V)
        zeta_w = 0.5 * (d + n)
        w0 = math.sqrt(k + d * n)
        wd = math.sqrt(w0 * w0 - zeta_w * zeta_w)
        A, B = a0, (-n * a0 + zeta_w * a0) / wd
        t, err, amax = 0.0, 0.0, 0.0
        num, ana = [], []
        while t < t_end - 1e-12:
            s = t5._rk4(t, s, dt, veh, aero, env)
            t += dt
            v = s[3:6]
            vn = math.sqrt(sum(x * x for x in v))
            cosang = sum(s[6 + i] * v[i] for i in range(3)) / vn
            sgn = 1.0 if (s[8] * v[1] - s[7] * v[2]) >= 0 else -1.0
            al = sgn * math.acos(max(-1.0, min(1.0, cosang)))
            an = math.exp(-zeta_w * t) * (A * math.cos(wd * t) + B * math.sin(wd * t))
            err = max(err, abs(al - an))
            amax = max(amax, abs(an))
            num.append(al)
            ana.append(an)
    finally:
        earth.gravity, earth.OMEGA = g_saved, om_saved
    return {"alpha0_deg": alpha0_deg, "mach": mach, "alt_km": h_km, "dt": dt,
            "omega_d_rad_s": wd, "period_s": 2 * math.pi / wd, "decay_rate_1_s": zeta_w,
            "max_abs_error_deg": math.degrees(err), "max_error_rel": err / a0}


if __name__ == "__main__":
    out = {"small": run(1.0), "large": run(5.0), "small_dt5": run(1.0, dt=5e-3)}
    with open(os.path.join(ROOT, "results", "weathercock_check.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
