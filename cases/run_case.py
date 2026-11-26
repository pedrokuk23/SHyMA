"""Runs paper test cases from the command line and checks them against the
published values: the same code and inputs as the reviewer workbench.

    python3 cases/run_case.py                      # all quick cases
    python3 cases/run_case.py --all                # also the slow ones (optimizer, full Monte Carlo)
    python3 cases/run_case.py cases/nominal_D.json # selected files
"""
import argparse
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.environ.setdefault("SHYMA_ROOT", ROOT)
sys.path.insert(0, os.path.join(ROOT, "web"))


def get(obj, path):
    for k in path.split("."):
        obj = obj[k]
    return obj


def run(case):
    import webdriver as wd
    fn = getattr(wd, case["function"])
    t0 = time.time()
    out = json.loads(fn(json.dumps(case["params"])) if case.get("params") is not None else fn())
    if case["function"] == "checks":
        out = {"items": out, "all_pass": all(x["pass"] for x in out)}
    rows, ok = [], True
    for e in case["expect"]:
        val = get(out, e["path"])
        if "max" in e:
            good = val <= e["max"]
            ref = "<= %g" % e["max"]
        elif isinstance(e["value"], bool):
            good, ref = val == e["value"], str(e["value"])
        else:
            good = abs(val - e["value"]) <= e["tol"]
            ref = "%.10g (tol %g)" % (e["value"], e["tol"])
        ok &= good
        rows.append((e["path"], val, ref, good))
    return ok, rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--all", action="store_true", help="include slow cases")
    a = ap.parse_args()
    files = a.files or sorted(glob.glob(os.path.join(HERE, "*.json")))
    n_fail = 0
    for f in files:
        case = json.load(open(f))
        if case.get("slow") and not (a.all or a.files):
            continue
        ok, rows, secs = run(case)
        n_fail += not ok
        print("%s  %s  [%s, %.0f s]" % ("PASS" if ok else "FAIL", case["name"], case["paper_reference"], secs))
        for path, val, ref, good in rows:
            print("    %-36s %-24s paper %s %s" % (path, val if not isinstance(val, float) else "%.10g" % val, ref,
                                                  "" if good else "<-- differs"))
    print("%d case(s) failed" % n_fail if n_fail else "all cases reproduce the paper")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
