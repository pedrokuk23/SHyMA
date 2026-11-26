"""Statistical tests of the held-out VSB-30 apogee residuals.

1) Calibration term. The held-out residuals share a common component induced
   by the uncertainty of (Isp, t_b1) calibrated on V01: their covariance is
   diag(sigma_i^2) + J_i C J_j^T, with J_i the apogee sensitivity of flight i
   to (Isp, t_b1) (central differences, booster calibrated on V01) and C the
   parameter covariance of run_vsb30_validation.py. The Mahalanobis statistic
   r^T S^-1 r is chi-square with n degrees of freedom.
2) Two-sided tests of the normalised residuals z_i = r_i / sigma_i against
   N(0, 1): upper tail (residuals larger than the predicted dispersion) and
   lower tail (smaller), and a t-test of the mean (bias).
3) Power of the upper-tail chi-square test for a model error sigma_d added in
   quadrature to the flight dispersion, and the one-sided upper 95 % bound of
   sigma_d (profile likelihood with r_i ~ N(0, sigma_i^2 + sigma_d^2)), with
   and without the calibration term.
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2, t as tdist

sys.path.insert(0, os.path.dirname(__file__))
import run_vsb30_validation as vv  # noqa: E402
from shyma import vehicle  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")


def per_flight_jacobian(rows, d_isp=1.0, d_b1=0.2):
    isp0, b10 = vehicle.ISP_VSB30, vehicle.S31_BURN
    J = []
    for r in rows:
        ap = vv.predict(r, isp0 + d_isp, b10)["apogee_km"]
        am = vv.predict(r, isp0 - d_isp, b10)["apogee_km"]
        bp = vv.predict(r, isp0, b10 + d_b1)["apogee_km"]
        bm = vv.predict(r, isp0, b10 - d_b1)["apogee_km"]
        J.append([(ap - am) / (2 * d_isp), (bp - bm) / (2 * d_b1)])
    return np.array(J)


def loglik(r, s2, sd):
    v = s2 + sd ** 2
    return -0.5 * float(np.sum(np.log(v) + r ** 2 / v))


def sigma_d_upper(r, s2, conf=0.95):
    """One-sided upper bound of sigma_d (km) from the profile likelihood (chi2_1 at 2*conf-1)."""
    grid = np.linspace(0.0, 60.0, 6001)
    ll = np.array([loglik(r, s2, g) for g in grid])
    best = ll.max()
    crit = chi2.ppf(2 * conf - 1, 1) / 2.0
    ok = grid[ll >= best - crit]
    return float(grid[int(np.argmax(ll))]), float(ok.max())


if __name__ == "__main__":
    v = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    cu = v["calibration_uncertainty"]
    Jv = np.array(cu["jacobian_km_per_s"])
    Ji = np.linalg.inv(Jv)
    sa, sr, rho = cu["obs_sigma_apogee_km"], cu["obs_sigma_range_km"], cu.get("obs_corr", 0.0)
    C = Ji @ np.array([[sa ** 2, rho * sa * sr], [rho * sa * sr, sr ** 2]]) @ Ji.T
    bl = [f for f in v["flights"] if not f["calibration"]]
    rows = {r["flight"]: r for r in vv.load()}
    Jf = per_flight_jacobian([rows[f["flight"]] for f in bl])
    r = np.array([f["apogee_model"] - f["apogee_actual"] for f in bl])
    s = np.array([f["mc_sigma_apogee_km"] for f in bl])
    n = len(bl)
    Ccal = Jf @ C @ Jf.T
    S = np.diag(s ** 2) + Ccal
    m2 = float(r @ np.linalg.solve(S, r))
    z = r / s
    q = float(np.sum(z ** 2))
    tt = float(np.mean(z) / (np.std(z, ddof=1) / np.sqrt(n)))
    sd_hat, sd_up = sigma_d_upper(r, s ** 2)
    # upper bound with the calibration term: profile over sigma_d with full covariance
    grid = np.linspace(0.0, 60.0, 3001)
    llc = []
    for g in grid:
        Sg = S + np.eye(n) * g ** 2
        sign, logdet = np.linalg.slogdet(Sg)
        llc.append(-0.5 * (logdet + float(r @ np.linalg.solve(Sg, r))))
    llc = np.array(llc)
    okc = grid[llc >= llc.max() - chi2.ppf(0.90, 1) / 2.0]
    apo_mean = float(np.mean([f["apogee_actual"] for f in bl]))
    crit = chi2.ppf(0.95, n)
    power = {}
    for frac in (0.5, 1.0, 1.5):           # sigma_d as a fraction of the mean flight sigma
        lam_scale = 1.0 + frac ** 2         # var ratio: (sigma^2 + sigma_d^2) / sigma^2 with sigma_d = frac*sigma
        power["%.1f" % frac] = float(chi2.sf(crit / lam_scale, n))
    out = {"n": n,
           "chi2_apogee_with_calibration": m2, "p": float(chi2.sf(m2, n)),
           "sigma_cal_apogee_km_per_flight": [float(np.sqrt(Ccal[i, i])) for i in range(n)],
           "sigma_cal_apogee_km": float(np.sqrt(np.mean(np.diag(Ccal)))),
           "jacobian_per_flight": Jf.tolist(),
           "chi2_z": q, "p_upper": float(chi2.sf(q, n)), "p_lower": float(chi2.cdf(q, n)),
           "p_two_sided": float(min(1.0, 2 * min(chi2.sf(q, n), chi2.cdf(q, n)))),
           "z_mean": float(np.mean(z)), "t_mean": tt, "p_mean": float(2 * tdist.sf(abs(tt), n - 1)),
           "power_upper_chi2_by_sigma_d_fraction": power,
           "sigma_d_hat_km": sd_hat, "sigma_d_upper95_km": sd_up,
           "sigma_d_upper95_pct": 100 * sd_up / apo_mean,
           "sigma_d_upper95_with_cal_km": float(okc.max()),
           "sigma_d_upper95_with_cal_pct": 100 * float(okc.max()) / apo_mean,
           "mean_flight_sigma_pct": float(np.mean(100 * s / np.array([f["apogee_actual"] for f in bl])))}
    for site in ("CLA", "Esrange"):
        idx = [i for i, f in enumerate(bl) if f["site"] == site]
        rs, Ss = r[idx], S[np.ix_(idx, idx)]
        qq = float(rs @ np.linalg.solve(Ss, rs))
        out["chi2_" + site] = qq
        out["p_" + site] = float(chi2.sf(qq, len(idx)))
    json.dump(out, open(os.path.join(ROOT, "results", "vsb30_zstats.json"), "w"), indent=1)
    print({k: x for k, x in out.items() if not isinstance(x, list)})
