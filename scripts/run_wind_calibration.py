"""Calibration of the residual-dispersion (residual-wind) amplitude k_w.

Target. Garcia et al. (JATM, 2011) state that the VSB-30 circular 3-sigma
dispersion radius is "around 50 km"; it is the design value obtained by
root-sum-square combination of the vehicle error sources, and all ten flown
impacts fell inside it (most within 2 sigma). It is therefore a design
dispersion consistent with flight experience, not an observed scatter.

Method. For the V01 configuration (CLA, 392 kg, 84 deg, 45 deg) and all
other dispersions of shyma.dispersion active, N dispersed flights are run
on a grid of k_w values with common random numbers (seed SEED). For each k_w
the 3-sigma radius (98.9 % containment) is estimated from the fitted
bivariate-normal impact covariance (low-variance estimator) and, for
comparison, as the empirical order statistic. k_w is found by linear
interpolation of the Gaussian estimate at 50 km, and its sampling interval
by bootstrap over the N draws. Because thrust and fin misalignment are not
modelled, k_w is a lumped residual-dispersion parameter expressed as an
equivalent residual wind.

Also: wind-weighting demonstration (8 m/s from 60 deg in the lowest 3 km).
"""
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import analysis, dispersion as D, earth, vehicle, wind  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
N = int(os.environ.get("MC_N_CAL", "600"))
SEED = 7
CFG = (392.0, 84.0, 45.0, "CLA")
TARGET_R = 50.0
K_GRID = [0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75]


def run_grid():
    rng = random.Random(SEED)
    base = [D.sample(rng, 0.0) for _ in range(N)]
    res = {}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as p:
        for k in K_GRID:
            draws = [dict(d, k_wind=k) for d in base]
            res[k] = p.map(D.one_flight, [(CFG, d) for d in draws])
            print("k=%.2f done" % k, flush=True)
    return res


def stats(out, idx=None):
    idx = range(len(out)) if idx is None else idx
    o = [out[i] for i in idx]
    n = len(o)
    lats, lons = [x["impact_lat"] for x in o], [x["impact_lon"] for x in o]
    e = D.ellipse(sum(lats) / n, sum(lons) / n, lats, lons)
    rm = sum(x["range_km"] for x in o) / n
    am = sum(x["apogee_km"] for x in o) / n
    return {"r3_gauss": e["r3sig_gauss_km"], "r3_emp": e["r3sig_empirical_km"],
            "sig_major": e["sigma_major_km"], "sig_minor": e["sigma_minor_km"],
            "range_sd": math.sqrt(sum((x["range_km"] - rm) ** 2 for x in o) / (n - 1)),
            "apogee_sd": math.sqrt(sum((x["apogee_km"] - am) ** 2 for x in o) / (n - 1))}


def solve_k(curve):
    ks = sorted(curve)
    for a, b in zip(ks[:-1], ks[1:]):
        ra, rb = curve[a], curve[b]
        if (ra - TARGET_R) * (rb - TARGET_R) <= 0 and ra != rb:
            return a + (TARGET_R - ra) * (b - a) / (rb - ra)
    return float("nan")


if __name__ == "__main__":
    t0 = time.time()
    res = run_grid()
    table = {k: stats(res[k]) for k in K_GRID}
    k_g = solve_k({k: table[k]["r3_gauss"] for k in K_GRID})
    k_e = solve_k({k: table[k]["r3_emp"] for k in K_GRID})
    # bootstrap interval of k_w (Gaussian estimator), resampling the draws
    r = random.Random(11)
    kb = []
    for _ in range(200):
        idx = [r.randrange(N) for _ in range(N)]
        kb.append(solve_k({k: stats(res[k], idx)["r3_gauss"] for k in K_GRID if k > 0}))
    kb = sorted(x for x in kb if x == x)
    # statistics at the calibrated k (interpolated between grid points)
    def at(key, k):
        ks = sorted(K_GRID)
        for a, b in zip(ks[:-1], ks[1:]):
            if a <= k <= b:
                return table[a][key] + (k - a) * (table[b][key] - table[a][key]) / (b - a)
        return float("nan")
    # wind-weighting demonstration
    mk = lambda: vehicle.vsb30(392.0)  # noqa: E731
    site = earth.STATIONS["CLA"]
    w8 = wind.uniform_wind(8.0, 60.0, top_m=3000.0)
    la0, lo0, tr0 = wind.impact_point(mk, site, 84.0, 45.0, None)
    la1, lo1, tr1 = wind.impact_point(mk, site, 84.0, 45.0, w8)
    ww = wind.wind_weighting(mk, site, 84.0, 45.0, w8)
    # unit-wind effect: impact shift for a 1 m/s wind from 60 deg, 0-3 km
    w1 = wind.uniform_wind(1.0, 60.0, top_m=3000.0)
    la2, lo2, _ = wind.impact_point(mk, site, 84.0, 45.0, w1)
    out = {"k_wind": k_g, "k_wind_ci95": [kb[int(0.025 * len(kb))], kb[int(0.975 * len(kb)) - 1]],
           "k_wind_empirical_estimator": k_e, "seed": SEED, "N": N, "k_grid": K_GRID,
           "target_r3sigma_km": TARGET_R,
           "table": {str(k): v for k, v in table.items()},
           "at_k": {key: at(key, k_g) for key in ("sig_major", "sig_minor", "range_sd", "apogee_sd",
                                                    "r3_emp", "r3_gauss")},
           "no_wind": table[0.0],
           "unit_wind_effect_km_per_mps": earth.great_circle_km(la0, lo0, la2, lo2),
           "weighting_demo": {"wind": "8 m/s from 60 deg, 0-3 km",
                              "unweighted_shift_km": earth.great_circle_km(la0, lo0, la1, lo1),
                              "unweighted_range_km": analysis.downrange_km(tr1, site)[-1],
                              "nominal_range_km": analysis.downrange_km(tr0, site)[-1],
                              "weighted_elevation": ww["elevation"], "weighted_azimuth": ww["azimuth"],
                              "residual_miss_m": ww["miss_m"], "iterations": len(ww["iterations"])},
           "wall_s": time.time() - t0}
    with open(os.path.join(ROOT, "results", "wind_calibration.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "table"}, indent=1))
