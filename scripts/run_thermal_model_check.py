"""Effect of the gas model on the aerothermal indicators: thermally perfect
air (vibrational excitation of N2 and O2, used in SHyMA) versus a
calorically perfect gas (gamma = 1.4), at the peak of the window for the
nominal concepts A-D (2 mm cylindrical leading edge)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
from shyma import analysis  # noqa: E402
from shyma.atmosphere import atmosphere  # noqa: E402

res = {}
for c in "ABCD":
    v, tr = cs.fly(c, 400.0, cs.NOMINAL[c], rec=1)
    w = cs.windows(c, v, tr)[0]
    ks = [k for k in range(len(tr.t)) if w["t_start"] <= tr.t[k] <= w["t_end"]]
    out = {"tp": [0.0, 0.0, 0.0], "cpg": [0.0, 0.0, 0.0]}
    for k in ks:
        T = atmosphere(tr.alt[k] / 1000.0)[0]
        qcw = analysis.sutton_graves(tr.rho[k], tr.v[k], cs.RN_LE, "cylinder")
        h0 = analysis.h_air(T) + 0.5 * tr.v[k] ** 2
        for key, fn, T0 in (("tp", analysis.wall_equilibrium, analysis.T_from_h(h0)),
                            ("cpg", analysis.wall_equilibrium_cpg, T + tr.v[k] ** 2 / (2 * 1004.5))):
            Tw, q = fn(qcw, T, tr.v[k])
            out[key] = [max(out[key][0], T0), max(out[key][1], Tw), max(out[key][2], q)]
    res[c] = {"T0_tp": out["tp"][0], "Tw_tp": out["tp"][1], "qhw_tp": out["tp"][2],
              "T0_cpg": out["cpg"][0], "Tw_cpg": out["cpg"][1], "qhw_cpg": out["cpg"][2]}
    print(c, res[c], flush=True)
json.dump(res, open(os.path.join(os.path.dirname(__file__), "..", "results", "thermal_model_check.json"), "w"),
          indent=1)
