"""Probabilistic threshold mass ceilings of the nominal settings.

The mass ceilings of run_optimizer.py and run_mass_budget.py are deterministic
(nominal trajectory, settings re-optimized at each mass). Here, for each concept
at its nominal setting and for the fixed-beta vehicle family, two outcomes are
recorded on a 20 kg mass grid with common random numbers (the same N dispersed
draws at every mass): the threshold window, and tau_use >= 3 s with a closed CLA
link (as in run_robust_selection.py).

A pointwise Clopper-Pearson bound at each grid mass is conservative (with
N = 200 it needs an observed probability of about 0.94). The probability is
therefore estimated between grid masses by linear interpolation of the observed
grid proportions (a saturated model, which needs no functional form), and its
pointwise one-sided 95 % lower band is the 5th percentile of the same
interpolation over a bootstrap of the draws, which keeps the pairing of the
common random numbers. The probabilistic ceiling is the largest mass
m >= 400 kg such that the lower band is >= P_MIN on [400 kg, m]; when the band
is below P_MIN at 400 kg, the range of masses in which it is >= P_MIN (if any)
is reported instead. A smooth logistic-quadratic fit, logit P = a + b x + c x^2,
is kept for comparison, with a Pearson goodness-of-fit statistic over the grid
(approximate: the grid masses share the draws). The outcome matrices are saved
in results/mass_ceiling_mc_outcomes.json.
"""
import json
import os
import random
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
import run_montecarlo as mc  # noqa: E402
from shyma import dispersion as D  # noqa: E402

N = int(os.environ.get("MCM_N", "200"))
N_BOOT = int(os.environ.get("MCM_BOOT", "300"))
SEED = 20261110
P_MIN = 0.90
TAU_REQ = 3.0
GRID = {"A": [300.0 + 20 * i for i in range(8)], "B": [380.0 + 20 * i for i in range(8)],
        "C": [380.0 + 20 * i for i in range(8)], "D": [320.0 + 20 * i for i in range(9)]}
METRICS = ("window", "use")
M_REF, M_SC = 400.0, 100.0


def one(args):
    c, m, d = args
    kw = {"coast": cs.COAST_D + d["dcoast"]} if c == "D" else {}
    v, tr = cs.fly(c, m * d["pl"], cs.NOMINAL[c] + d["del"], az=cs.AZIMUTH + d["daz"], rho_factor=d["rho"],
                   dt=0.02, rec=5, wind=D.wind_of(d), dT=D.temp_of(d), beta_tv=cs.BETA_TV * d["beta"], beta_stack=cs.BETA_TV,
                   cd_scale=d["cd"], impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                   burn_factor=(d["burn1"], d["burn2"]), **kw)
    w = cs.windows(c, v, tr, max(cs.DSEP + d["dsep"], 0.0))[0]
    tu = cs.useful_window(tr, w, ant_loss_db=mc.ANT_LOSS_DB, fade_db=mc.FADE_DB, mask_deg=mc.ELEV_MASK_DEG)
    return int(bool(w["met"])), int(tu >= TAU_REQ)


def fit(masses, Y, idx):
    """ML logistic-quadratic fit; Y[j][i] outcome of draw i at mass j; idx: draws used."""
    x = np.repeat((np.array(masses) - M_REF) / M_SC, len(idx))
    y = np.concatenate([np.array(Y[j])[idx] for j in range(len(masses))]).astype(float)

    def nll(p):
        eta = p[0] + p[1] * x + p[2] * x * x
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta)) + 1e-6 * float(p @ p)
    return minimize(nll, np.zeros(3), method="BFGS").x


def curve(p, mm):
    x = (mm - M_REF) / M_SC
    return 1.0 / (1.0 + np.exp(-(p[0] + p[1] * x + p[2] * x * x)))


def interp_curve(masses, Y, idx, mm):
    """Observed proportions at the grid masses (draws idx), linearly interpolated."""
    pg = np.array([np.array(Y[j])[idx].mean() for j in range(len(masses))])
    return np.interp(mm, np.array(masses), pg)


def pearson_gof(masses, Y, p):
    """Pearson chi^2 of the quadratic fit over the grid and its approximate p-value."""
    from scipy.stats import chi2
    n = len(Y[0])
    e = curve(p, np.array(masses))
    o = np.array([np.mean(Y[j]) for j in range(len(masses))])
    x2 = float(np.sum(n * (o - e) ** 2 / np.maximum(e * (1 - e), 1e-9)))
    dof = len(masses) - 3
    return x2, dof, float(chi2.sf(x2, dof))


def ceiling(mm, band):
    ok = band >= P_MIN
    i0 = int(np.argmin(np.abs(mm - M_REF)))
    if ok[i0]:
        j = i0
        while j + 1 < len(mm) and ok[j + 1]:
            j += 1
        return {"ceiling_kg": float(mm[j]), "band_kg": None}
    if ok.any():
        sel = mm[ok]
        return {"ceiling_kg": None, "band_kg": [float(sel.min()), float(sel.max())]}
    return {"ceiling_kg": None, "band_kg": None}


if __name__ == "__main__":
    t0 = time.time()
    rng = random.Random(SEED)
    k = mc.k_wind_calibrated()
    ds = [D.sample(rng, k) for _ in range(N)]
    out = {"n": N, "n_boot": N_BOOT, "seed": SEED, "p_min": P_MIN, "model": "linear interpolation of grid proportions; bootstrap pointwise lower band",
           "concepts": {}}
    rb = np.random.default_rng(7)
    outcomes = {}
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as pool:
        for c, ms in GRID.items():
            Y = {key: [] for key in METRICS}
            rows = []
            for m in ms:
                res = pool.map(one, [(c, m, d) for d in ds], chunksize=4)
                r = {"mass_kg": m}
                for key, col in zip(METRICS, zip(*res)):
                    Y[key].append(list(col))
                    met = sum(col)
                    lo, hi = D.clopper_pearson(met, N)
                    r[key] = {"p": met / N, "p_lo95": lo, "p_hi95": hi}
                rows.append(r)
                print(c, m, r["window"]["p"], r["use"]["p"], flush=True)
            mm = np.arange(ms[0], ms[-1] + 0.5, 1.0)
            rec = {"grid": rows}
            for key in METRICS:
                allidx = np.arange(N)
                boot_idx = [rb.integers(0, N, N) for _ in range(N_BOOT)]
                est = interp_curve(ms, Y[key], allidx, mm)
                band = np.percentile(np.array([interp_curve(ms, Y[key], b, mm) for b in boot_idx]), 5.0, axis=0)
                p_q = fit(ms, Y[key], allidx)
                band_q = np.percentile(np.array([curve(fit(ms, Y[key], b), mm) for b in boot_idx]), 5.0, axis=0)
                x2, dof, pgof = pearson_gof(ms, Y[key], p_q)
                i400 = int(np.argmin(np.abs(mm - 400)))
                rec[key] = {"p_at_400": float(est[i400]), "band_at_400": float(band[i400]),
                            "band_max": float(band.max()), "mass_at_band_max": float(mm[int(np.argmax(band))]),
                            **ceiling(mm, band),
                            "curve_masses": mm[::10].tolist(), "curve_p": est[::10].tolist(), "curve_band": band[::10].tolist(),
                            "quadratic": {"coef": p_q.tolist(), "p_at_400": float(curve(p_q, mm)[i400]),
                                          "band_at_400": float(band_q[i400]), **ceiling(mm, band_q),
                                          "pearson_chi2": x2, "dof": dof, "p_gof": pgof}}
                print(c, key, {kk: rec[key][kk] for kk in ("ceiling_kg", "band_kg", "band_at_400", "band_max")},
                      "quadratic gof p=%.3f" % pgof, flush=True)
            outcomes[c] = {"masses": ms, **{key: Y[key] for key in METRICS}}
            out["concepts"][c] = rec
    out["cpu_s"] = time.time() - t0
    json.dump(out, open(os.path.join(mc.ROOT, "results", "mass_ceiling_mc.json"), "w"), indent=1)
    json.dump(outcomes, open(os.path.join(mc.ROOT, "results", "mass_ceiling_mc_outcomes.json"), "w"))
