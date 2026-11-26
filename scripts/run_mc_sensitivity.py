"""Which dispersions drive the dwell time? Spearman rank correlations between
each dispersed input and the dwell of the 1000-run Monte Carlo set of each
concept (results/montecarlo.json, "with_wind").

The draws are regenerated from the seed with the same random stream as
run_montecarlo.run(), so no trajectory is recomputed; one run per concept is
re-flown to check that the regenerated draws are the ones that were flown.
The residual-wind profile enters through its two unit vectors (low and upper
level), projected on the downrange (azimuth 45 deg) and crossrange axes.
The squared rank correlations, normalised to sum to one, rank the inputs by
association with the dwell (a screening measure; not a variance
decomposition and not Sobol indices).
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_montecarlo as mc  # noqa: E402
from shyma import dispersion as D  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
LABELS = {
    "imp1": "S31 impulse", "imp2": "S30 impulse", "burn1": "S31 burn time", "burn2": "S30 burn time",
    "dry": "stage dry mass", "pl": "payload mass", "cd": "axial drag", "rho": "pressure scale",
    "dT0": "temperature offset, low", "dT1": "temperature offset, high",
    "del": "launcher elevation", "daz": "launcher azimuth", "beta": "test-vehicle ballistic coeff.",
    "dsep": "separation delay (A)", "dcoast": "coast timer (D)",
    "w0d": "residual wind, low, downrange", "w0c": "residual wind, low, crossrange",
    "w1d": "residual wind, upper, downrange", "w1c": "residual wind, upper, crossrange",
}


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = 0.5 * (i + j)
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    sy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return sxy / (sx * sy) if sx > 0 and sy > 0 else 0.0


def wind_components(d):
    r = random.Random(d["wseed"])
    u0 = (r.gauss(0, 1), r.gauss(0, 1))
    u1 = (r.gauss(0, 1), r.gauss(0, 1))
    c = math.sqrt(0.5)
    return {"w0d": c * (u0[0] + u0[1]), "w0c": c * (u0[0] - u0[1]),
            "w1d": c * (u1[0] + u1[1]), "w1c": c * (u1[0] - u1[1])}


def main():
    res = json.load(open(os.path.join(ROOT, "results", "montecarlo.json")))
    n, k = res["N"], res["k_wind"]
    out = {"n": n, "concepts": {}}
    for c in "ABCD":
        rng = random.Random(mc.SEED)   # same draws for every concept (run_montecarlo.run)
        draws = [D.sample(rng, k, 1.0) for _ in range(n)]
        y = res["with_wind"][c]["window_samples"]
        assert len(y) == n
        check = mc.one((c, draws[0]))["window_s"]
        assert abs(check - y[0]) < 1e-9, (c, check, y[0])
        cols = {key: [d[key] for d in draws] for key in LABELS if key in draws[0]}
        wc = [wind_components(d) for d in draws]
        for key in ("w0d", "w0c", "w1d", "w1c"):
            cols[key] = [w[key] for w in wc]
        rho = {key: spearman(v, y) for key, v in cols.items()}
        if c != "A":
            rho.pop("dsep", None)
        if c != "D":
            rho.pop("dcoast", None)
        tot = sum(r * r for r in rho.values())
        share = {key: r * r / tot for key, r in rho.items()}
        top = sorted(rho, key=lambda key: -abs(rho[key]))
        out["concepts"][c] = {"rho": rho, "share": share, "ranking": top}
        print(c, [(LABELS[key], round(rho[key], 2), round(share[key], 2)) for key in top[:5]], flush=True)
    out["labels"] = LABELS
    with open(os.path.join(ROOT, "results", "mc_sensitivity.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
