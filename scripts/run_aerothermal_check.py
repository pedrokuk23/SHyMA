"""Code-to-code check of the stagnation heating: Sutton-Graves (used by
SHyMA) against the Fay-Riddell correlation for a sphere, at the window-entry
and window-exit states of the nominal concepts.

Fay-Riddell (equilibrium boundary layer, negligible dissociation since
T0 < 2500 K, so the Lewis-number term vanishes):
  q = 0.763 Pr^-0.6 (rho_w mu_w / rho_e mu_e)^0.1 sqrt(rho_e mu_e du/dx) (h_0 - h_w),
with the Newtonian velocity gradient du/dx = (1/R_n) sqrt(2 (p_e - p_inf)/rho_e),
p_e the pitot pressure (Rayleigh, gamma = 1.4), h_0 and T_0 from thermally
perfect air, Sutherland viscosity and Pr = 0.71. The wall is cold (300 K)
for the cold-wall comparison and at the SHyMA radiative-equilibrium
temperature for the hot-wall comparison.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import analysis  # noqa: E402
from shyma.atmosphere import atmosphere, R_AIR  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
G = 1.4
PR = 0.71


def mu(T):
    return 1.458e-6 * T ** 1.5 / (T + 110.4)


def pitot(p, M):
    a = ((G + 1) ** 2 * M * M / (4 * G * M * M - 2 * (G - 1))) ** (G / (G - 1))
    return p * a * (1 - G + 2 * G * M * M) / (G + 1)


def fay_riddell(alt_km, M, rn, Tw):
    T, p, rho, a = atmosphere(alt_km)
    V = M * a
    h0 = analysis.h_air(T) + 0.5 * V * V
    T0 = analysis.T_from_h(h0)
    pe = pitot(p, M)
    rhoe, mue = pe / (R_AIR * T0), mu(T0)
    rhow, muw = pe / (R_AIR * Tw), mu(Tw)
    dudx = math.sqrt(2 * (pe - p) / rhoe) / rn
    hw = analysis.h_air(Tw)
    return 0.763 * PR ** -0.6 * (rhow * muw / (rhoe * mue)) ** 0.1 * math.sqrt(rhoe * mue * dudx) * (h0 - hw), h0, hw


def sutton_graves_hw(alt_km, M, rn, Tw):
    T, p, rho, a = atmosphere(alt_km)
    V = M * a
    q = analysis.sutton_graves(rho, V, rn, "sphere")
    h0 = analysis.h_air(T) + 0.5 * V * V
    return q, q * (1 - analysis.h_air(Tw) / h0)


if __name__ == "__main__":
    cs = json.load(open(os.path.join(ROOT, "results", "casestudy.json")))
    out = []
    for c, n in cs["nominal"].items():
        for i, tag in ((0, "entry"), (1, "exit")):
            h, M = n["window_alt_km"][i], n["window_mach"][i]
            for rn in (0.002, 0.05):
                qsg_cw, _ = sutton_graves_hw(h, M, rn, 300.0)
                qfr_cw, h0, hw = fay_riddell(h, M, rn, 300.0)
                T, p, rho, a = atmosphere(h)
                tw, qsg_hw = analysis.wall_equilibrium(qsg_cw, T, M * a, 0.85)
                qfr_hw, _, _ = fay_riddell(h, M, rn, tw)
                out.append({"concept": c, "state": tag, "alt_km": h, "mach": M, "rn_m": rn,
                            "sg_cold_MW": qsg_cw * (1 - hw / h0) / 1e6, "fr_cold_MW": qfr_cw / 1e6,
                            "sg_hot_MW": qsg_hw / 1e6, "fr_hot_MW": qfr_hw / 1e6, "Tw_K": tw,
                            "ratio_cold": qsg_cw * (1 - hw / h0) / qfr_cw, "ratio_hot": qsg_hw / qfr_hw})
                print(out[-1], flush=True)
    rc = [o["ratio_cold"] for o in out]
    rh = [o["ratio_hot"] for o in out]
    res = {"cases": out, "ratio_cold_range": [min(rc), max(rc)], "ratio_hot_range": [min(rh), max(rh)]}
    with open(os.path.join(ROOT, "results", "aerothermal_check.json"), "w") as f:
        json.dump(res, f, indent=1)
    print(res["ratio_cold_range"], res["ratio_hot_range"])
