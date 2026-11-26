"""Transient lumped-capacitance estimate of the leading-edge temperature.

The radiative-equilibrium wall temperature of the main text is an upper bound.
Here a solid cylindrical leading edge of radius R (2 mm, Inconel 718,
rho = 8190 kg/m^3, c = 435 J/kg/K) is treated as a lumped mass per unit span
heated on the stagnation line (f = 1, upper bound) or with a circumferential
average of the hot-wall flux (f = 0.3) and radiating over its perimeter:

    rho c R / 2  dT/dt = f q_hw(T) - eps sigma T^4,  q_hw = q_cw (1 - h(T)/h0)

from the exposure time (shroud jettison / separation) of each nominal concept.
Conduction into the structure is neglected (it lowers T further) and so is the
latent heat: the model is not valid above the solidus, whose crossing time is
reported. The output gives the temperature at the end of the window, its ratio
to the equilibrium value, the thermal time constant at the window peak flux
(with dh/dT of thermally perfect air) and the heat loads from exposure to the
window end: cold-wall, transient-wall and radiative-equilibrium hot-wall.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
from shyma import analysis  # noqa: E402
from shyma.atmosphere import atmosphere  # noqa: E402

RHO, CP, EPS, T_INIT = 8190.0, 435.0, 0.85, 300.0
T_SOLIDUS = 1533.0


def cp_inconel(T):
    return CP + 0.20 * (T - 300.0)


def dhdt_air(T):
    return (analysis.h_air(T + 5.0) - analysis.h_air(T - 5.0)) / 10.0
FACTORS = (1.0, 0.3)


def run(c):
    v, tr = cs.fly(c, 400.0, cs.NOMINAL[c], rec=1)
    w = cs.windows(c, v, tr)[0]
    if w["t_start"] is None:
        return None
    if c == "A":
        t_exp = v.t_burnout + cs.DSEP
    elif c == "D":
        t_exp = w["t_start"] - 2.0
    else:
        t_exp = 60.0
    cap0 = RHO * cs.RN_LE / 2.0
    out = {"t_expose": t_exp, "window": [w["t_start"], w["t_end"]]}
    for fct in FACTORS:
        T, Tw_end, Teq_end, Tmax, Teq_max = T_INIT, None, None, T_INIT, 0.0
        tau_peak, qpk = None, -1.0
        Qcw = Qtr = Qhw = 0.0
        t_sol = None
        for k in range(1, len(tr.t)):
            t0, t1 = tr.t[k - 1], tr.t[k]
            if t1 <= t_exp or t0 > w["t_end"]:
                continue
            a = tr.alt[k]
            qc = analysis.sutton_graves(tr.rho[k], tr.v[k], cs.RN_LE, "cylinder") if a < 120e3 else 0.0
            Tinf = atmosphere(max(a, 0.0) / 1000.0)[0]
            h0 = analysis.h_air(Tinf) + 0.5 * tr.v[k] ** 2
            dt = t1 - max(t0, t_exp)
            nsub = max(1, int(dt / 0.002))
            for _ in range(nsub):
                qh = fct * qc * max(1.0 - analysis.h_air(T) / h0, 0.0)
                Qtr += qh * dt / nsub
                T += (qh - EPS * analysis.SIGMA * T ** 4) / (cap0 * cp_inconel(T)) * dt / nsub
            if t_sol is None and T >= T_SOLIDUS:
                t_sol = t1 - w["t_start"]
            teq, qeq = analysis.wall_equilibrium(fct * qc, Tinf, tr.v[k], EPS)
            Qcw += fct * qc * dt
            Qhw += qeq * dt
            if w["t_start"] <= t1 <= w["t_end"]:
                Tmax, Teq_max = max(Tmax, T), max(Teq_max, teq)
                if qc > qpk:
                    qpk = qc
                    # linearized time constant around equilibrium: cap / (4 eps sigma Teq^3 + f q_cw cp/h0)
                    tau_peak = cap0 * cp_inconel(teq) / (4 * EPS * analysis.SIGMA * teq ** 3 + fct * qc * dhdt_air(teq) / h0)
                Tw_end, Teq_end = T, teq
        out["f%.1f" % fct] = {"T_end_K": Tw_end, "Teq_end_K": Teq_end, "T_max_window_K": Tmax,
                              "Teq_max_window_K": Teq_max, "tau_s": tau_peak,
                              "ratio_rise": (Tmax - T_INIT) / max(Teq_max - T_INIT, 1.0),
                              "t_solidus_from_window_start_s": t_sol,
                              "Q_cold_MJ": Qcw / 1e6, "Q_transient_MJ": Qtr / 1e6, "Q_hot_eq_MJ": Qhw / 1e6}
    return out


if __name__ == "__main__":
    res = {"material": {"rho": RHO, "cp_300K": CP, "cp_slope": 0.20, "eps": EPS, "R_m": cs.RN_LE,
                        "T_init": T_INIT, "T_solidus": T_SOLIDUS}, "concepts": {}}
    for c in "ABCD":
        res["concepts"][c] = run(c)
        print(c, res["concepts"][c], flush=True)
    json.dump(res, open(os.path.join(os.path.dirname(__file__), "..", "results", "le_transient.json"), "w"), indent=1)
