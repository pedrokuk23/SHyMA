"""Validation against the ten VSB-30 flights of Garcia et al. (JATM, 2011).

1) Held-out prediction: model calibrated on V01 (shyma/calibration.json),
   flights V02-V10 predicted from payload mass, site, elevation and azimuth
   only (WGS-84 Earth; the spherical Earth of the earlier ASMAT satellite tool is also
   reported, with its own V01 calibration).
2) Statistics with uncertainty: bootstrap interval of RMS(model) -
   RMS(operators), Wilcoxon signed-rank test on |errors|, per-site results.
3) Leave-one-out: the model is recalibrated on each flight in turn and the
   other nine are predicted (distribution of held-out RMS errors).
4) Calibration uncertainty: Jacobian of (apogee, range) w.r.t. (Isp, t_b1)
   and the parameter covariance induced by the flight-to-flight dispersion
   predicted by the Monte Carlo for the V01 configuration.
5) Normalised residuals: for every flight configuration the Monte Carlo
   (all dispersions, calibrated k_w) gives sigma_apogee and sigma_range;
   z = (model - actual)/sigma is tested against N(0,1).
The operators' pre-flight nominal predictions of the same table are used as
a reference; note that they were made before flight with nominal motor
data, while SHyMA benefits from a post-flight calibration on V01.
"""
import csv
import json
import math
import os
import random
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import analysis, calibration, dispersion as D, earth, flight, vehicle  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
N_MC = int(os.environ.get("MC_N_VAL", "300"))


def load():
    rows = []
    with open(os.path.join(ROOT, "data", "vsb30_flights_garcia2011.csv")) as f:
        for r in csv.DictReader(line for line in f if not line.startswith("#")):
            rows.append(r)
    return rows


def cfg(r):
    return {"flight": r["flight"], "payload_kg": float(r["payload_kg"]),
            "elevation": float(r["elevation_deg"]), "azimuth": float(r["azimuth_deg"]),
            "site": r["site"], "apogee_km": float(r["apogee_actual_km"]),
            "range_km": float(r["range_actual_km"])}


def predict(r, isp=None, b1=None):
    v = vehicle.vsb30(float(r["payload_kg"]), isp=isp, s31_burn=b1)
    site = earth.STATIONS[r["site"]]
    tr = flight.fly(v, site, float(r["elevation_deg"]), float(r["azimuth_deg"]), 60.0,
                    record_every=1)
    ta, ha, _ = tr.apogee()
    t100 = 0.0
    for k in range(1, len(tr.t)):
        a0, a1 = tr.alt[k - 1], tr.alt[k]
        if a0 > 100e3 and a1 > 100e3:
            t100 += tr.t[k] - tr.t[k - 1]
        elif (a0 - 100e3) * (a1 - 100e3) < 0:
            w = (100e3 - a0) / (a1 - a0)
            t100 += (1 - w) * (tr.t[k] - tr.t[k - 1]) if a1 > 100e3 else w * (tr.t[k] - tr.t[k - 1])
    return {"apogee_km": ha / 1e3, "range_km": analysis.downrange_km(tr, site)[-1],
            "t_above_100km": t100}


def stats(e):
    n = len(e)
    m = sum(e) / n
    rms = math.sqrt(sum(x * x for x in e) / n)
    sd = math.sqrt(sum((x - m) ** 2 for x in e) / (n - 1)) if n > 1 else float("nan")
    return {"mean": m, "rms": rms, "sd": sd, "max_abs": max(abs(x) for x in e), "n": n}


def held_out_set(rows, isp=None, b1=None, cal="V01"):
    recs = []
    for r in rows:
        p = predict(r, isp, b1)
        a_act, R_act = float(r["apogee_actual_km"]), float(r["range_actual_km"])
        rec = dict(flight=r["flight"], site=r["site"], payload_kg=float(r["payload_kg"]),
                   elevation=float(r["elevation_deg"]), azimuth=float(r["azimuth_deg"]),
                   apogee_actual=a_act, apogee_nominal=float(r["apogee_nominal_km"]),
                   apogee_model=p["apogee_km"], range_actual=R_act,
                   range_nominal=float(r["range_nominal_km"]), range_model=p["range_km"],
                   t100_actual=float(r["time_above_100km_s"]), t100_model=p["t_above_100km"],
                   calibration=(r["flight"] == cal))
        rec["apogee_err_model_pct"] = 100 * (rec["apogee_model"] - a_act) / a_act
        rec["apogee_err_nominal_pct"] = 100 * (rec["apogee_nominal"] - a_act) / a_act
        rec["range_err_model_km"] = rec["range_model"] - R_act
        rec["range_err_nominal_km"] = rec["range_nominal"] - R_act
        rec["t100_err_model_pct"] = 100 * (rec["t100_model"] - rec["t100_actual"]) / rec["t100_actual"]
        recs.append(rec)
    return recs


def summary(recs):
    held = [f for f in recs if not f["calibration"]]
    out = {"apogee_model_pct": stats([f["apogee_err_model_pct"] for f in held]),
           "apogee_nominal_pct": stats([f["apogee_err_nominal_pct"] for f in held]),
           "range_model_km": stats([f["range_err_model_km"] for f in held]),
           "range_nominal_km": stats([f["range_err_nominal_km"] for f in held]),
           "t100_model_pct": stats([f["t100_err_model_pct"] for f in held])}
    for site in ("CLA", "Esrange"):
        b = [f for f in held if f["site"] == site]
        if b:
            out["apogee_model_pct_" + site] = stats([f["apogee_err_model_pct"] for f in b])
    return out


def paired_tests(recs):
    from scipy.stats import wilcoxon
    held = [f for f in recs if not f["calibration"]]
    em = [f["apogee_err_model_pct"] for f in held]
    en = [f["apogee_err_nominal_pct"] for f in held]
    r = random.Random(3)
    n = len(em)
    diffs = []
    for _ in range(20000):
        idx = [r.randrange(n) for _ in range(n)]
        rm = math.sqrt(sum(em[i] ** 2 for i in idx) / n)
        rn = math.sqrt(sum(en[i] ** 2 for i in idx) / n)
        diffs.append(rm - rn)
    diffs.sort()
    w = wilcoxon([abs(a) for a in em], [abs(b) for b in en])
    return {"rms_diff_pct": stats(em)["rms"] - stats(en)["rms"],
            "rms_diff_ci95": [diffs[int(0.025 * len(diffs))], diffs[int(0.975 * len(diffs)) - 1]],
            "p_model_worse": sum(d > 0 for d in diffs) / len(diffs),
            "wilcoxon_abs_p": float(w.pvalue)}


def mc_config(args):
    (payload, el, az, site), seed, k = args
    rng = random.Random(seed)
    out = [D.one_flight(((payload, el, az, site), D.sample(rng, k))) for _ in range(N_MC // 2)]
    return out


def loo(rows):
    res = []
    for r in rows:
        f = cfg(r)
        c = calibration.calibrate(f)
        recs = held_out_set(rows, c["isp"], c["b1"], cal=r["flight"])
        s = summary(recs)
        res.append({"calibrated_on": r["flight"], "isp": c["isp"], "b1": c["b1"],
                    "apogee_rms_pct": s["apogee_model_pct"]["rms"],
                    "apogee_mean_pct": s["apogee_model_pct"]["mean"],
                    "range_rms_km": s["range_model_km"]["rms"]})
        print("LOO", res[-1], flush=True)
    return res


if __name__ == "__main__":
    rows = load()
    out = {"isp": vehicle.ISP_VSB30, "s31_burn": vehicle.S31_BURN, "earth_model": earth.MODEL}
    out["flights"] = held_out_set(rows)
    out["stats_held_out"] = summary(out["flights"])
    out["tests"] = paired_tests(out["flights"])
    print(json.dumps(out["stats_held_out"], indent=1), out["tests"], flush=True)

    # spherical Earth (geodesy of the earlier ASMAT satellite tool), recalibrated on V01
    earth.set_model("sphere")
    cs = calibration.calibrate(calibration.V01)
    sph = held_out_set(rows, cs["isp"], cs["b1"])
    out["sphere"] = {"isp": cs["isp"], "b1": cs["b1"], "stats_held_out": summary(sph)}
    earth.reset_model()

    # leave-one-out
    out["loo"] = loo(rows)

    # calibration uncertainty and normalised residuals (needs k_w)
    with open(os.path.join(ROOT, "results", "wind_calibration.json")) as f:
        wc = json.load(f)
    k = wc["k_wind"]
    J = json.load(open(os.path.join(ROOT, "shyma", "calibration.json")))["jacobian_km_per_s"]
    tasks = []
    for i, r in enumerate(rows):
        c = (float(r["payload_kg"]), float(r["elevation_deg"]), float(r["azimuth_deg"]), r["site"])
        tasks += [(c, 1000 + 2 * i, k), (c, 1001 + 2 * i, k)]
    with Pool(int(os.environ.get("MC_PROCS", "2"))) as p:
        mc = p.map(mc_config, tasks)
    zs = []
    for i, f in enumerate(out["flights"]):
        o = mc[2 * i] + mc[2 * i + 1]
        n = len(o)
        am = sum(x["apogee_km"] for x in o) / n
        rm = sum(x["range_km"] for x in o) / n
        sa = math.sqrt(sum((x["apogee_km"] - am) ** 2 for x in o) / (n - 1))
        sr = math.sqrt(sum((x["range_km"] - rm) ** 2 for x in o) / (n - 1))
        f["mc_sigma_apogee_km"], f["mc_sigma_range_km"] = sa, sr
        f["z_apogee"] = (f["apogee_model"] - f["apogee_actual"]) / sa
        f["z_range"] = (f["range_model"] - f["range_actual"]) / sr
        if not f["calibration"]:
            zs.append((f["z_apogee"], f["z_range"]))
        else:   # V01: observation covariance of the calibration (apogee, range), with correlation
            cov = sum((x["apogee_km"] - am) * (x["range_km"] - rm) for x in o) / (n - 1)
            obs = (sa, sr, cov / (sa * sr))
    # calibration uncertainty: the dispersion of the calibration flight V01 (apogee and range,
    # with their correlation, all dispersions of the Monte Carlo model including the
    # residual-dispersion profile, which V01 itself experienced) mapped through the Jacobian
    C = calibration.parameter_covariance(J, *obs)
    C0 = calibration.parameter_covariance(J, wc["at_k"]["apogee_sd"], wc["at_k"]["range_sd"])
    C1 = calibration.parameter_covariance(J, wc["no_wind"]["apogee_sd"], wc["no_wind"]["range_sd"])
    out["calibration_uncertainty"] = {"jacobian_km_per_s": J, "sigma_isp_s": math.sqrt(C[0][0]),
                                      "sigma_b1_s": math.sqrt(C[1][1]),
                                      "corr": C[0][1] / math.sqrt(C[0][0] * C[1][1]),
                                      "obs_sigma_apogee_km": obs[0], "obs_sigma_range_km": obs[1],
                                      "obs_corr": obs[2],
                                      "diagonal_wind_calibration": {"sigma_isp_s": math.sqrt(C0[0][0]),
                                                                    "sigma_b1_s": math.sqrt(C0[1][1])},
                                      "without_residual_profile": {"sigma_isp_s": math.sqrt(C1[0][0]),
                                                                   "sigma_b1_s": math.sqrt(C1[1][1])}}
    from scipy.stats import chi2
    za = sum(z[0] ** 2 for z in zs)
    zr = sum(z[1] ** 2 for z in zs)
    nb = len(zs)
    out["normalised_residuals"] = {
        "n": nb, "chi2_apogee": za, "p_apogee": float(chi2.sf(za, nb)),
        "chi2_range": zr, "p_range": float(chi2.sf(zr, nb)),
        "rms_z_apogee": math.sqrt(za / nb), "rms_z_range": math.sqrt(zr / nb),
        "mean_z_apogee": sum(z[0] for z in zs) / nb, "mean_z_range": sum(z[1] for z in zs) / nb,
        "k_wind": k, "n_mc_per_flight": N_MC}
    print(json.dumps(out["normalised_residuals"], indent=1))
    with open(os.path.join(ROOT, "results", "vsb30_validation.json"), "w") as f:
        json.dump(out, f, indent=2)
