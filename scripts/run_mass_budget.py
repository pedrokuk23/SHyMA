"""Mass margin of each architecture at its nominal launch settings, for two
families of separated test vehicles.

  fixed-beta     : the ballistic coefficient stays at BETA_TV whatever the
                   mass (a family of vehicles scaled with the mass; the
                   assumption of the paper's mass sweeps);
  fixed-geometry : the shape and drag curve are those of the 400 kg vehicle,
                   so beta = BETA_TV * m / 400 (mass added to the same vehicle,
                   e.g. batteries, mechanisms or instrumentation).
Concepts C and D keep the test vehicle on the stack, so both families
coincide for them.

For each family the interval of test-vehicle masses around 400 kg (1 kg
resolution) is found for which the nominal trajectory still gives (a) the threshold window and (b) a
useful test time (window and CLA link, run_casestudy.useful_window) of
at least TAU_REQ. The margin above 400 kg is the mass that services,
instrumentation and growth allowance may take in each architecture before the
requirement is lost, at fixed launch settings (re-tuning the settings raises
it, see the mass ceilings of run_optimizer.py).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
TAU_REQ = 3.0
M_REF = 400.0


def metrics(c, m, family):
    beta = cs.BETA_TV * (m / M_REF if family == "fixed-geometry" else 1.0)
    kw = {"coast": cs.COAST_D} if c == "D" else {}
    v, tr = cs.fly(c, m, cs.NOMINAL[c], beta_tv=beta, dt=0.02, rec=2, **kw)
    w = cs.windows(c, v, tr)[0]
    return {"window_s": w["duration"], "met": bool(w["met"]), "tau_use": cs.useful_window(tr, w)}


def bounds(c, family, key):
    """Feasible mass interval around the 400 kg reference: bisection upwards
    for the ceiling and downwards (to 250 kg) for the floor."""
    ok = (lambda r: r["met"]) if key == "window" else (lambda r: r["met"] and r["tau_use"] >= TAU_REQ)
    if not ok(metrics(c, M_REF, family)):
        return None, None

    def bisect(a, b):          # a feasible, b infeasible
        while abs(b - a) > 1.0:
            m = 0.5 * (a + b)
            if ok(metrics(c, m, family)):
                a = m
            else:
                b = m
        return a
    hi = 700.0 if ok(metrics(c, 700.0, family)) else bisect(M_REF, 700.0)
    lo = 250.0 if ok(metrics(c, 250.0, family)) else bisect(M_REF, 250.0)
    return lo, hi


def main():
    out = {"tau_req": TAU_REQ, "m_ref": M_REF, "beta_ref": cs.BETA_TV, "concepts": {}}
    for c in "ABCD":
        fams = ("fixed-beta", "fixed-geometry") if c in "AB" else ("fixed-beta",)
        rec = {}
        for fam in fams:
            r = {"at_400": metrics(c, M_REF, fam),
                 "sweep": {str(m): metrics(c, float(m), fam) for m in range(300, 520, 20)}}
            for key in ("window", "useful"):
                lo, hi = bounds(c, fam, key)
                r["floor_" + key + "_kg"], r["ceiling_" + key + "_kg"] = lo, hi
                r["margin_" + key + "_kg"] = None if hi is None else hi - M_REF
            rec[fam] = r
            print(c, fam, r["floor_window_kg"], r["ceiling_window_kg"], r["floor_useful_kg"],
                  r["ceiling_useful_kg"], flush=True)
        out["concepts"][c] = rec
    with open(os.path.join(ROOT, "results", "mass_budget.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
