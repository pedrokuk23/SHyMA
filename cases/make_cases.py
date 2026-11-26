"""Writes the paper test cases (cases/*.json) from the published results.

Each case holds the inputs of one run of the workbench functions
(web/webdriver.py: nominal, optimize, montecarlo, vsb30, nesc, checks) and the
values that run must reproduce. The same files can be run from the command line
(python3 cases/run_case.py) or loaded in the reviewer workbench
(web/shyma_workbench.html, "Load a paper case").
"""
import json
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)
import run_casestudy as cs          # noqa: E402
import run_montecarlo as mc         # noqa: E402
from shyma import dispersion as D   # noqa: E402


def load(name):
    with open(os.path.join(ROOT, "results", name + ".json")) as f:
        return json.load(f)


WINDOW = {"mach": list(cs.WINDOW["mach"]), "alt_km": list(cs.WINDOW["alt_km"]), "q_kpa": list(cs.WINDOW["q_kpa"]),
          "min_duration": cs.WINDOW["min_duration"]}
OBJECTIVE = {"mach_min": cs.OBJECTIVE["mach"][0], "q_min": cs.OBJECTIVE["q_kpa"][0]}


def base(concept, el):
    return {"concept": concept, "payload": 400.0, "el": el, "az": cs.AZIMUTH, "beta": cs.BETA_TV, "rn_le": cs.RN_LE,
            "dsep": cs.DSEP, "coast": cs.COAST_D, "window": WINDOW, "objective": OBJECTIVE, "phase": "default"}


def write(name, case):
    with open(os.path.join(ROOT, "cases", name + ".json"), "w") as f:
        json.dump(case, f, indent=1)
    print("cases/%s.json" % name)


def main():
    nom = load("casestudy")["nominal"]
    for c in "ABCD":
        n = nom[c]
        write("nominal_%s" % c, {
            "name": "Nominal concept %s (400 kg, paper inputs)" % c, "paper_reference": "Table 4, Sec. 6.2",
            "function": "nominal", "tab": "nominal", "runtime": "about 5 s", "slow": False,
            "params": base(c, cs.NOMINAL[c]),
            "expect": [{"path": "summary." + k, "value": n[k], "tol": 1e-6}
                       for k in ("window_s", "objective_s", "steady_s", "apogee_km", "impact_range_km") if n.get(k) is not None]})
    opt = load("optimizer")
    for c in "ABCD":
        b = opt["max_dwell_400"][c]
        free = {"el": None, "coast": None, "dsep": None, "az": None}
        free.update(opt["free"][c])
        p = base(c, 68.0)
        p.update(target="threshold", objective_fn="max_dwell", free=free, constraints={})
        exp = [{"path": "best.threshold", "value": b["threshold_s"], "tol": 1e-6}]
        exp += [{"path": "best.x." + k, "value": v, "tol": 1e-9} for k, v in b["x"].items()]
        write("optimizer_%s" % c, {
            "name": "Inverse design, concept %s: longest threshold dwell at 400 kg" % c, "paper_reference": "Table 5",
            "function": "optimize", "tab": "opt", "runtime": "2 to 4 min", "slow": True, "params": p, "expect": exp})
    mcj = load("montecarlo")
    for c in "ABCD":
        p = {"concept": c, "n": 20, "seed": mc.SEED, "wind": True, "scale": 1.0,
             "sigma": dict(D.SIG)}
        write("montecarlo_%s_first20" % c, {
            "name": "Monte Carlo, concept %s: first 20 runs of the paper set" % c, "paper_reference": "Table 6, Fig. S9",
            "function": "montecarlo", "tab": "mc", "runtime": "about 1 min", "slow": False, "params": p,
            "expect": [{"path": "comparison.max_abs_diff_s", "value": 0.0, "tol": 1e-9}]})
        w = mcj["with_wind"][c]
        p = dict(p, n=mcj["N"])
        write("montecarlo_%s_full" % c, {
            "name": "Monte Carlo, concept %s: full 1000-run paper set" % c, "paper_reference": "Table 6",
            "function": "montecarlo", "tab": "mc", "runtime": "20 to 60 min", "slow": True, "params": p,
            "expect": [{"path": "comparison.max_abs_diff_s", "value": 0.0, "tol": 1e-9},
                       {"path": "result.p_window", "value": w["p_window"], "tol": 1e-12},
                       {"path": "result.window_p05", "value": w["window_p05"], "tol": 1e-9},
                       {"path": "result.window_p50", "value": w["window_p50"], "tol": 1e-9}]})
    v = load("vsb30_validation")["stats_held_out"]
    write("vsb30_validation", {
        "name": "VSB-30: ten flights, V01 calibration, nine held-out", "paper_reference": "Fig. 2, Table S19",
        "function": "vsb30", "tab": "val", "runtime": "about 15 s", "slow": False, "params": None,
        "expect": [{"path": "summary.apogee_model_pct." + k, "value": v["apogee_model_pct"][k], "tol": 1e-6}
                   for k in ("mean", "rms", "max_abs")] +
                  [{"path": "summary.range_model_km.rms", "value": v["range_model_km"]["rms"], "tol": 1e-6}]})
    write("formula_checks", {
        "name": "Formula checks against independent references", "paper_reference": "Secs. 3-4, Supplementary S1",
        "function": "checks", "tab": "checks", "runtime": "under 1 s", "slow": False, "params": None,
        "expect": [{"path": "all_pass", "value": True}]})
    ne = load("nesc_checkcases")
    write("nesc_checkcases", {
        "name": "NASA NESC atmospheric check-cases 5, 9, 10 and 17", "paper_reference": "Table S2, Fig. S1",
        "function": "nesc", "tab": "val", "runtime": "about 5 s", "slow": False, "params": None,
        "note": "The workbench compares with a down-sampled copy of the NASA reference data; the bounds below are the "
                "paper's acceptance criterion (maximum error below the NASA spread for cases 5, 9 and 10).",
        "expect": [{"path": "%s.alt.max_abs_err" % k, "max": ne[k]["metrics"]["alt_m"]["max_spread"]}
                   for k in ("case05", "case09", "case10")]})


if __name__ == "__main__":
    main()
