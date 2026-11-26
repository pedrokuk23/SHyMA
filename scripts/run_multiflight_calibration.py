"""Multi-flight calibration of the VSB-30 booster with a launch-site atmosphere.

The paper calibrates (Isp_vac, t_b1) on one flight (V01). Here both parameters
are fitted by weighted least squares (Gauss-Newton, central differences) to the
apogees and ranges of all ten flights, each flown with NRLMSISE-00 for its site
and date (12 UTC, F10.7 = F10.7a = 100, Ap = 5), with the flight-to-flight
standard deviations of the Monte Carlo model (results/vsb30_validation.json) as
weights. The parameter covariance is (J^T W J)^-1, inflated by the reduced
chi-square when it exceeds one. Leave-one-out refits give the cross-validated
apogee error, and the nominal concepts are flown with the fitted parameters
(with US76 and with the December NRLMSISE-00 atmosphere of CLA).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
import run_montecarlo as mc  # noqa: E402
import run_vsb30_validation as vv  # noqa: E402
from shyma import analysis, earth, flight, vehicle  # noqa: E402
from shyma.atmosphere import msis_atmosphere  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
_ATM = {}


def atm_of(r):
    key = r["flight"]
    if key not in _ATM:
        la, lo = earth.STATIONS[r["site"]][:2]
        _ATM[key] = msis_atmosphere(la, lo, r["date"] + "T12:00", f107=100.0, f107a=100.0, ap=5.0)
    return _ATM[key]


def predict(r, isp, b1):
    v = vehicle.vsb30(float(r["payload_kg"]), isp=isp, s31_burn=b1)
    site = earth.STATIONS[r["site"]]
    tr = flight.fly(v, site, float(r["elevation_deg"]), float(r["azimuth_deg"]), 60.0, record_every=5, atm=atm_of(r))
    return np.array([tr.apogee()[1] / 1e3, analysis.downrange_km(tr, site)[-1]])


def residuals(rows, sig, x):
    pred = [predict(r, x[0], x[1]) for r in rows]
    obs = [np.array([float(r["apogee_actual_km"]), float(r["range_actual_km"])]) for r in rows]
    return np.concatenate([(o - p) / s for o, p, s in zip(obs, pred, sig)]), pred


def fit(rows, sig, x0, iters=6):
    x = np.array(x0, float)
    for _ in range(iters):
        e0, _ = residuals(rows, sig, x)
        J = np.zeros((len(e0), 2))
        for j, h in enumerate((1.0, 0.2)):
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            J[:, j] = -(residuals(rows, sig, xp)[0] - residuals(rows, sig, xm)[0]) / (2 * h)
        dx = np.linalg.lstsq(J, e0, rcond=None)[0]
        x = x + dx
        if abs(dx[0]) < 0.01 and abs(dx[1]) < 0.002:
            break
    e, pred = residuals(rows, sig, x)
    return x, J, e, pred


if __name__ == "__main__":
    t0 = time.time()
    rows = vv.load()
    vj = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    sg = {f["flight"]: (f["mc_sigma_apogee_km"], f["mc_sigma_range_km"]) for f in vj["flights"]}
    sig = [np.array(sg[r["flight"]]) for r in rows]
    x0 = (vehicle.ISP_VSB30, vehicle.S31_BURN)
    x, J, e, pred = fit(rows, sig, x0)
    chi2 = float(e @ e)
    dof = len(e) - 2
    C = np.linalg.inv(J.T @ J) * max(1.0, chi2 / dof)
    apo_err = [100 * (p[0] - float(r["apogee_actual_km"])) / float(r["apogee_actual_km"]) for r, p in zip(rows, pred)]
    loo = []
    for i in range(len(rows)):
        rr = rows[:i] + rows[i + 1:]
        ss = sig[:i] + sig[i + 1:]
        xi, _, _, _ = fit(rr, ss, x, iters=4)
        p = predict(rows[i], xi[0], xi[1])
        loo.append({"flight": rows[i]["flight"], "isp": float(xi[0]), "b1": float(xi[1]),
                    "apogee_err_pct": 100 * (p[0] - float(rows[i]["apogee_actual_km"])) / float(rows[i]["apogee_actual_km"])})
        print(loo[-1], flush=True)
    # nominal concepts with the fitted parameters
    conc = {}
    for label, atm in (("us76", None), ("msis_cla", mc.msis_cla())):
        conc[label] = {}
        for c in "ABCD":
            v, tr = cs.fly(c, 400.0, cs.NOMINAL[c], atm=atm, isp=float(x[0]), s31_burn=float(x[1]))
            s = cs.summarize(c, v, tr, light=True, atm=atm)
            conc[label][c] = {"window_s": s["window_s"], "burnout_mach": s["burnout_mach"]}
    loo_err = [l["apogee_err_pct"] for l in loo]
    out = {"isp": float(x[0]), "b1": float(x[1]), "sigma_isp": float(np.sqrt(C[0, 0])), "sigma_b1": float(np.sqrt(C[1, 1])),
           "corr": float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])), "chi2": chi2, "dof": dof,
           "apogee_err_pct": apo_err, "apogee_rms_pct": float(np.sqrt(np.mean(np.square(apo_err)))),
           "loo": loo, "loo_apogee_rms_pct": float(np.sqrt(np.mean(np.square(loo_err)))),
           "loo_sigma_isp": float(np.std([l["isp"] for l in loo], ddof=1)),
           "loo_sigma_b1": float(np.std([l["b1"] for l in loo], ddof=1)),
           "v01_isp": x0[0], "v01_b1": x0[1], "concepts": conc, "wall_s": time.time() - t0}
    print({k: v for k, v in out.items() if k not in ("loo",)})
    json.dump(out, open(os.path.join(ROOT, "results", "multiflight_calibration.json"), "w"), indent=1)
