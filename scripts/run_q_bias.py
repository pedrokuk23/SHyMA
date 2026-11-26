"""Effect of the HyShot 2 dynamic-pressure bias on the re-entry concepts: the
dwell of Concepts B and C with the density scaled by 1/1.08 and 1/1.19 (the
over-prediction of q at 29 and 35 km in the HyShot 2 reconstruction)."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import run_casestudy as cs  # noqa: E402

if __name__ == "__main__":
    out = {"factors": {}, "max_abs_change_s": 0.0}
    base = {}
    for fac in (1.0, 1 / 1.08, 1 / 1.19):
        row = {}
        for c in "BC":
            v, tr = cs.fly(c, 400.0, cs.NOMINAL[c], rho_factor=fac)
            row[c] = cs.summarize(c, v, tr, light=True)["window_s"]
        out["factors"]["%.4f" % fac] = row
        if fac == 1.0:
            base = row
        else:
            out["max_abs_change_s"] = max(out["max_abs_change_s"], *(abs(row[c] - base[c]) for c in "BC"))
    json.dump(out, open(os.path.join(os.path.dirname(__file__), "..", "results", "q_bias.json"), "w"), indent=1)
    print(out)
