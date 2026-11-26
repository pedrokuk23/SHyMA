"""Builds web/shyma_workbench.html: a single-file reviewer workbench that runs
the SHyMA code itself in the browser (Pyodide) and compares the results with
the published ones. Usage: python3 web/build_html.py"""
import base64
import csv
import hashlib
import io
import json
import os
import subprocess
import zipfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R = lambda *p: os.path.join(ROOT, *p)  # noqa: E731


def load(name):
    with open(R("results", name)) as f:
        return json.load(f)


def published():
    cs = load("casestudy.json")
    nom = {c: {k: v for k, v in d.items() if k != "profile"} for c, d in cs["nominal"].items()}
    mcj = load("montecarlo.json")
    mc = {}
    for c, d in mcj["with_wind"].items():
        e = dict(d)
        e["xy"] = [[round(x, 3) for x in d["xy"][0]], [round(y, 3) for y in d["xy"][1]]]
        e.pop("failures", None)
        mc[c] = e
    meta = {k: mcj[k] for k in ("N", "sigma", "k_wind", "seed", "ant_loss_db", "fade_db", "nominal_elevation")}
    meta["no_wind"] = {c: {k: d[k] for k in ("n", "p_window", "window_p05", "window_p50")} for c, d in mcj["no_wind"].items()}
    meta["doubled_sigma"] = {c: {k: d[k] for k in ("n", "p_window", "window_p05", "window_p50")}
                             for c, d in mcj["doubled_sigma"].items()}
    v = load("vsb30_validation.json")
    ne = load("nesc_checkcases.json")
    nesc = {k: {"metrics": ne[k]["metrics"], "ref_profiles": ne[k]["ref_profiles"]} for k in ne}
    with open(R("validation_data", "Atmos_17_TwoStageRocketToOrbit", "Atmos_17_sim_05.csv")) as f:
        rows = list(csv.DictReader(f))
    r0 = min(rows, key=lambda r: abs(float(r["time"]) - 20.0))
    keep = ["time"] + ["eiPosition_ft_" + k for k in "XYZ"] + ["eiVelocity_ft_s_" + k for k in "XYZ"]
    nesc["case17"]["init_row"] = {k: r0[k] for k in keep}
    opt = load("optimizer.json")
    return {"optimizer": {"free": opt["free"], "max_dwell_400": opt["max_dwell_400"], "objective_350": opt["objective_350"],
                          "mass_ceiling": {c: {"ceiling_kg": v.get("ceiling_kg")} for c, v in opt["mass_ceiling"].items()}},
            "nominal": nom, "montecarlo": mc, "montecarlo_meta": meta, "k_wind": load("wind_calibration.json")["k_wind"],
            "vsb30": {"flights": v["flights"], "stats_held_out": v["stats_held_out"]}, "nesc": nesc}


FILES = (["shyma/" + f for f in sorted(os.listdir(R("shyma"))) if f.endswith((".py", ".json"))]
         + ["scripts/run_casestudy.py", "scripts/run_montecarlo.py", "scripts/run_vsb30_validation.py",
            "scripts/run_nesc_checkcases.py", "scripts/launch_optimizer.py", "data/vsb30_flights_garcia2011.csv", "results/wind_calibration.json",
            "web/webdriver.py"])


def main():
    pub = json.dumps(published(), separators=(",", ":"))
    buf = io.BytesIO()
    manifest = []
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in FILES:
            data = open(R(f), "rb").read()
            z.writestr(f, data)
            manifest.append({"file": f, "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
        z.writestr("web/published.json", pub)
    try:
        rev = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True).stdout.strip()
    except Exception:
        rev = ""
    code_hash = hashlib.sha256("".join(m["sha256"] for m in manifest).encode()).hexdigest()[:12]
    info = {"code_revision": rev or "v2.4", "code_hash": code_hash, "files": manifest}
    html = open(R("web", "template.html")).read()
    html = html.replace("__ZIP_B64__", base64.b64encode(buf.getvalue()).decode())
    html = html.replace("__BUILD_INFO__", json.dumps(info))
    cases = []
    for f in sorted(os.listdir(R("cases"))) if os.path.isdir(R("cases")) else []:
        if f.endswith(".json"):
            cases.append(json.load(open(R("cases", f))))
    order = {"nominal": 0, "optimize": 1, "montecarlo": 2, "vsb30": 3, "nesc": 4, "checks": 5}
    cases.sort(key=lambda c: (order.get(c["function"], 9), c.get("slow", False), c["name"]))
    html = html.replace("__CASES__", json.dumps(cases, separators=(",", ":")))
    out = R("web", "shyma_workbench.html")
    with open(out, "w") as f:
        f.write(html)
    print(out, round(os.path.getsize(out) / 1e6, 2), "MB")


if __name__ == "__main__":
    main()
