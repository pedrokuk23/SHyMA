"""Verification of the geometry-based stability derivatives (shyma.stability).

1. Exact limits of the component methods:
   * slender cone: C_N_alpha = 2 (base area) at 2/3 of the length;
   * Diederich/Barrowman panel slope: 2 pi for A -> inf, pi A / 2 for A -> 0;
   * supersonic panel slope: 4/beta for A -> inf; continuity at the sonic-leading-edge
     condition beta cot(Lambda) = 1;
   * K_W(B): 1 at tau = 0, 2 at tau = 1; K_W(B) + K_B(W) = (1 + tau)^2;
   * complete elliptic integral E(k) against reference values.
2. Inverse dynamic check (cf. Greenberg, NACA TN 2340; Jaffe, JPL TR 32-757): the
   pseudo-5-DOF integrator of shyma.trajectory5 is run for a coasting stack at
   constant conditions with the geometric derivatives of the VSB-30 at Mach 3; the
   damped oscillation of the angle of attack is fitted and inverted for C_m_alpha
   and C_mq, which must return the input values.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import drag as Dr, earth, geometry as G, stability as St, trajectory5 as t5, vehicle  # noqa: E402
from shyma.trajectory import Env  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def limits():
    out = {}
    cone = G.Configuration([(0.0, 0.0), (3.0, 0.5)], [], 1.0)
    s = St.static(cone, 2.0)
    out["cone_cna"] = s["cn_alpha"]
    out["cone_xcp_over_L"] = s["x_cp"] / 3.0
    big = G.FinSet("a", 0.0, 1.0, 1.0, 1000.0, 0.0)
    small = G.FinSet("b", 0.0, 1.0, 1.0, 0.001, 0.0)
    out["diederich_large_A_over_2pi"] = St.panel_cna_subsonic(big, 0.0) / (2 * math.pi)
    out["diederich_small_A_over_piA2"] = St.panel_cna_subsonic(small, 0.0) / (math.pi * small.aspect / 2.0)
    out["supersonic_large_A_times_beta_over_4"] = St.panel_cna_supersonic(big, 2.0) * math.sqrt(3.0) / 4.0
    f = G.FinSet("c", 0.0, 1.0, 0.0, 1.0, 1.0)            # 45 deg LE sweep: sonic LE at beta = 1
    m0 = math.sqrt(2.0)
    out["sonic_le_jump"] = St.panel_cna_supersonic(f, m0 + 1e-6) / St.panel_cna_supersonic(f, m0 - 1e-6) - 1.0
    out["kwb_0"], out["kwb_1"] = St.k_wb(1e-9), St.k_wb(0.998)
    out["kwb_plus_kbw_minus_1ptau2"] = max(abs(St.k_wb(t) + St.k_bw(t) - (1 + t) ** 2) for t in (0.1, 0.3, 0.5, 0.8))
    out["E_0"], out["E_05"], out["E_1"] = St.ellip_e(0.0), St.ellip_e(0.5), St.ellip_e(1.0)
    out["E_05_ref"] = 1.4674622093394272       # E(k = 0.5)
    return out


def inverse_check(mach=3.0, h_km=20.0, alpha0_deg=1.0, dt=1e-3, t_end=4.0):
    T, p, rho, a = 216.65, 5474.9, 0.088035, 295.07          # US76 at 20 km
    env = Env(atm=lambda alt_km, level, rf: (T, p, rho, a))
    veh = vehicle.vsb30(400.0, cd_scale=0.0)
    for s in veh.stages:                                     # coasting full stack
        s.t_ign, s.burn_time, s.t_drop = 1e9, 1.0, 2e9
    g = G.vsb30_geometry()
    items = G.vsb30_mass_items(veh, g["stations"])
    m, xcg, I = G.mass_properties(items, 0.0)
    d = St.pitch_derivatives(g["full"], mach, xcg)
    aero = t5.Aero5(0.557, cn_alpha=d["cn_alpha"], margin_cal=d["margin_cal"], cmq=d["cmq"])
    aero.inertia = lambda veh_, t, tc=None: I
    gs, om = earth.gravity, earth.OMEGA
    earth.gravity = lambda x, y, z: (0.0, 0.0, 0.0)
    earth.OMEGA = 0.0
    try:
        V = mach * a
        a0 = math.radians(alpha0_deg)
        s = [earth.RE + h_km * 1e3, 0.0, 0.0, 0.0, V, 0.0, 0.0, math.cos(a0), math.sin(a0), 0.0, 0.0, 0.0]
        t, ts, al = 0.0, [], []
        while t < t_end - 1e-12:
            s = t5._rk4(t, s, dt, veh, aero, env)
            t += dt
            v = s[3:6]
            vn = math.sqrt(sum(x * x for x in v))
            c = sum(s[6 + i] * v[i] for i in range(3)) / vn
            sg = 1.0 if (s[8] * v[1] - s[7] * v[2]) >= 0 else -1.0
            ts.append(t)
            al.append(sg * math.acos(max(-1.0, min(1.0, c))))
    finally:
        earth.gravity, earth.OMEGA = gs, om
    lam, w, A, B, rms = St.fit_damped_oscillation(ts, al)
    q = 0.5 * rho * V * V
    inv = St.derivatives_from_oscillation(lam, w, q, V, veh.ref_area, 0.557, m, I, d["cn_alpha"])
    return {"mach": mach, "input": {"cn_alpha": d["cn_alpha"], "cm_alpha": d["cm_alpha"], "cmq": d["cmq"]},
            "fit": {"decay_1_s": lam, "omega_rad_s": w, "rms_rad": rms},
            "recovered": inv,
            "err_cm_alpha_pct": 100 * (inv["cm_alpha"] / d["cm_alpha"] - 1),
            "err_cmq_pct": 100 * (inv["cmq"] / d["cmq"] - 1)}


def drag_checks():
    """Gas-dynamics and friction building blocks of shyma.drag against reference values."""
    out = {}
    out["pm_nu_M2_deg"], out["pm_nu_M2_ref"] = math.degrees(Dr.pm_nu(2.0)), 26.380   # NACA Rep. 1135
    b, M2, p21, _ = Dr.oblique(2.0, math.radians(10.0))
    out["oblique_M2_d10_beta_deg"], out["oblique_beta_ref"] = math.degrees(b), 39.31  # NACA Rep. 1135
    out["oblique_M2_d10_p21"], out["oblique_p21_ref"] = p21, 1.7066
    cp, Ms, bs = Dr.taylor_maccoll(2.0, math.radians(10.0))
    out["cone10_M2_shock_deg"], out["cone10_shock_ref"] = math.degrees(bs), 31.2         # NACA Rep. 1135 cone charts
    # hypersonic small-disturbance cone pressure (Rasmussen) as an independent check at Mach 6 and 8
    for M, th in ((6.0, 10.0), (8.0, 5.0)):
        t_ = math.radians(th)
        K = M * t_
        hsd = t_ * t_ * (1 + ((Dr.GAM + 1) * K * K + 2) / ((Dr.GAM - 1) * K * K + 2) * math.log((Dr.GAM + 1) / 2 + 1 / (K * K)))
        out["cone%d_M%d_cp_over_hsd" % (th, M)] = Dr.taylor_maccoll(M, t_)[0] / hsd
    out["cf_KS_Re1e7"], out["cf_KS_ref"] = Dr.cf_incompressible(1e7), 0.002934          # Karman-Schoenherr
    out["cf_VD_ratio_M2_Re1e7"] = Dr.cf_van_driest(1e7, 2.0) / Dr.cf_incompressible(1e7)
    return out


if __name__ == "__main__":
    out = {"limits": limits(), "inverse": inverse_check(), "drag": drag_checks()}
    with open(os.path.join(ROOT, "results", "stability_check.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
