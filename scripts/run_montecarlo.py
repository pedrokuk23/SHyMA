"""Monte Carlo dispersion analysis of the four nominal concepts (A-D).

Dispersions: shyma.dispersion (impulse, burn time, dry mass, payload mass,
drag, pressure/density scale, altitude-correlated temperature offset,
launcher elevation and azimuth, test-vehicle ballistic coefficient,
Concept A separation delay, and the calibrated residual-wind profile).
Runs: with residual wind (N per concept), without wind (N/2), and with all
1-sigma values doubled (N_2X), to expose the sensitivity to the assumed
sigmas; with a model-form impulse error (N/2); with the calibration and k_w
uncertainty, one sampled calibration per run (N/2); and with a 25 % test-vehicle
ballistic-coefficient dispersion for Concepts A and D (N/2); with the NRLMSISE-00
December atmosphere of CLA as the baseline (N/2); with a model discrepancy equal to
the upper 95 % bound of the held-out apogee model error (N/2); with the booster
parameters sampled from the leave-one-out spread and from the multi-flight
calibration (N/2 each; the latter, fitted with the site atmospheres, is
flown with the NRLMSISE-00 atmosphere of CLA); and, for Concept D, with a random attitude error of the
stack at S30 ignition (1 and 2 deg per axis, N/2 each). All sets use common random
numbers: every concept and every set restart the same random streams, so run i of
any set uses the aleatory draws of run i of the base set (paired differences, exact
McNemar test; the draws are also common to the four concepts).
The useful test time tau_use (continuous interval in the window with a closed
CLA link, run_casestudy.useful_window) is recorded for every run. Statistics carry 95 % confidence intervals (Clopper-Pearson for
probabilities, bootstrap for percentiles and ellipse axes).

Link closure: the S-band margin of the nominal budget (0 dBi vehicle
antenna) is degraded by ANT_LOSS_DB (aft-oblique aspect, antenna pattern)
and FADE_DB (multipath at low elevation); the link closes if the degraded
minimum CLA margin during the test is >= 0 dB and the CLA elevation stays above the 3 deg mask.
"""
import json
import math
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
from shyma import dispersion as D  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
N = int(os.environ.get("MC_N", "1000"))
N_2X = int(os.environ.get("MC_N_2X", "300"))
SEED = 20261006
ANT_LOSS_DB = 10.0
FADE_DB = 6.0
ELEV_MASK_DEG = 3.0   # visibility mask also required for link closure


_MSIS = None


def msis_cla():
    """NRLMSISE-00 profile of CLA for a December campaign (as in run_casestudy.py)."""
    global _MSIS
    if _MSIS is None:
        from shyma.atmosphere import msis_atmosphere
        _MSIS = msis_atmosphere(cs.SITE[0], cs.SITE[1], "2021-12-14T15:00", f107=100.0, f107a=100.0, ap=5.0)
    return _MSIS


def k_wind_calibrated():
    with open(os.path.join(ROOT, "results", "wind_calibration.json")) as f:
        return json.load(f)["k_wind"]


def one(args):
    concept, d = args
    kw = {"coast": cs.COAST_D + d["dcoast"]} if concept == "D" else {}
    if "isp" in d:          # calibration-uncertainty set: per-run booster parameters
        kw.update(isp=d["isp"], s31_burn=d["b1"])
    atm = msis_cla() if d.get("msis") else None
    if "att" in d:          # attitude error at S30 ignition (Concept D)
        kw["ign_att"] = d["att"]
    if "tip" in d:          # tip-off attitude error just after rail exit (run_tipoff_check.py)
        kw["tip_att"] = d["tip"]
    v, tr = cs.fly(concept, 400.0 * d["pl"], cs.NOMINAL[concept] + d["del"],
                   az=cs.AZIMUTH + d["daz"], rho_factor=d["rho"], dt=0.02, rec=5,
                   wind=D.wind_of(d), dT=D.temp_of(d), beta_tv=cs.BETA_TV * d["beta"], beta_stack=cs.BETA_TV,
                   cd_scale=d["cd"], impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                   burn_factor=(d["burn1"], d["burn2"]), atm=atm, **kw)
    dsep = max(cs.DSEP + d["dsep"], 0.0)
    s = cs.summarize(concept, v, tr, dsep=dsep, light=True, atm=atm)
    w = cs.windows(concept, v, tr, dsep)[0]
    tau_use = cs.useful_window(tr, w, ant_loss_db=ANT_LOSS_DB, fade_db=FADE_DB, mask_deg=ELEV_MASK_DEG)
    keep = ("window_s", "window_met", "objective_s", "objective_met", "steady_s", "apogee_km",
            "impact_range_km", "impact_lat", "impact_lon", "burnout_mach", "burnout_alt_km",
            "t_window", "qalpha_max_kpa_deg", "q_max_ascent_kpa")
    out = {k: s.get(k) for k in keep}
    out["tau_use"] = tau_use
    if "window_stations" in s:
        out["cla_el_min"] = s["window_stations"]["CLA"]["el_min"]
        out["cla_margin_min"] = s["window_stations"]["CLA"]["margin_min_db"]
        out["window_mach_min"] = min(s["window_mach"])
    return out


def analyse(out):
    n = len(out)
    lats = [o["impact_lat"] for o in out]
    lons = [o["impact_lon"] for o in out]
    e = D.ellipse(sum(lats) / n, sum(lons) / n, lats, lons)
    w = [o["window_s"] for o in out]
    met = sum(o["window_met"] for o in out)
    link = sum(1 for o in out if o["window_met"] and o.get("cla_margin_min") is not None
               and o.get("cla_el_min") is not None and o["cla_el_min"] >= ELEV_MASK_DEG
               and o["cla_margin_min"] - ANT_LOSS_DB - FADE_DB >= 0.0)
    obj = sum(o["objective_met"] for o in out)
    use = sum(o["tau_use"] >= 3.0 for o in out)
    apo = [o["apogee_km"] for o in out]
    am = sum(apo) / n
    els = [o["cla_el_min"] for o in out if o.get("cla_el_min") is not None]
    qa = [o["qalpha_max_kpa_deg"] for o in out if o.get("qalpha_max_kpa_deg") is not None]
    p = lambda q: D.percentile(w, q)  # noqa: E731
    res = {"n": n, "p_window": met / n, "p_window_ci95": D.clopper_pearson(met, n),
           "p_window_link": link / n, "p_window_link_ci95": D.clopper_pearson(link, n),
           "p_objective": obj / n, "p_objective_ci95": D.clopper_pearson(obj, n),
           "p_use": use / n, "p_use_ci95": D.clopper_pearson(use, n),
           "met": "".join("1" if o["window_met"] else "0" for o in out),
           "use_ok": "".join("1" if o["tau_use"] >= 3.0 else "0" for o in out),
           "tau_use_p05": D.percentile([o["tau_use"] for o in out], 0.05),
           "tau_use_p50": D.percentile([o["tau_use"] for o in out], 0.5),
           "window_p05": p(0.05), "window_p50": p(0.5), "window_p95": p(0.95),
           "window_p05_ci95": D.bootstrap_ci(w, lambda x: D.percentile(x, 0.05)),
           "steady_p50": D.percentile([o["steady_s"] for o in out], 0.5),
           "apogee_mean": am, "apogee_std": math.sqrt(sum((a - am) ** 2 for a in apo) / (n - 1)),
           "range_mean": sum(o["impact_range_km"] for o in out) / n,
           "ellipse": {k: v for k, v in e.items() if k not in ("x_km", "y_km")},
           "xy": [e["x_km"], e["y_km"]],
           "window_samples": w}
    xs, ys = e["x_km"], e["y_km"]
    res["ellipse"]["sigma_major_ci95"] = D.bootstrap_ci(
        list(zip(lats, lons)), lambda s: D.ellipse(sum(a for a, _ in s) / len(s), sum(b for _, b in s) / len(s),
                                                   [a for a, _ in s], [b for _, b in s])["sigma_major_km"],
        n_boot=200)
    if els:
        res["cla_el_min_p05"] = D.percentile(els, 0.05)
    if qa:
        res["qalpha_p99"] = D.percentile(qa, 0.99)
    fails = [o for o in out if not o["window_met"]]
    res["failures"] = [{k: o.get(k) for k in ("window_s", "burnout_mach", "window_mach_min", "t_window")}
                       for o in fails[:20]]
    _ = xs, ys
    return res


def model_form_sigma():
    """Common (both stages) impulse factor that represents the model-form error:
    conservatively, the whole held-out apogee RMS error of the VSB-30 validation is
    attributed to the model. The Monte Carlo applies the factor to the
    propellant mass at fixed Isp (imp1, imp2), so it is converted with the
    sensitivity of the V01 apogee to that same factor, computed here by central
    differences with the calibrated booster."""
    from shyma import calibration, earth, flight, vehicle
    v = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    rms = v["stats_held_out"]["apogee_model_pct"]["rms"] / 100.0
    f0 = calibration.V01
    site = earth.STATIONS[f0["site"]]

    def apo(fac):
        veh = vehicle.vsb30(f0["payload_kg"], impulse_factor=(fac, fac))
        tr = flight.fly(veh, site, f0["elevation"], f0["azimuth"], 60.0, record_every=5,
                        model=calibration.MODEL, dt=calibration.DT)
        return tr.apogee()[1] / 1000.0
    h = 0.01
    a0, ap, am = apo(1.0), apo(1.0 + h), apo(1.0 - h)
    sens = (ap - am) / (2.0 * h) / a0          # relative apogee change per relative impulse change
    return rms / sens


def calibration_uncertainty():
    """Joint 1-sigma uncertainty of the V01 calibration (Isp, t_b1, with
    correlation) and of the residual-dispersion amplitude k_w (from the
    interval of k_w obtained with a +-10 % uncertainty on the 50 km target)."""
    from shyma import vehicle
    v = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))["calibration_uncertainty"]
    w = json.load(open(os.path.join(ROOT, "results", "wind_calibration.json")))
    tab = {float(k): x for k, x in w["table"].items()}
    ks = sorted(tab)

    def k_at(target):
        for a0, b0 in zip(ks[:-1], ks[1:]):
            ra, rb = tab[a0]["r3_gauss"], tab[b0]["r3_gauss"]
            if (ra - target) * (rb - target) <= 0:
                return a0 + (target - ra) * (b0 - a0) / (rb - ra)
        raise ValueError(target)
    return {"isp": vehicle.ISP_VSB30, "b1": vehicle.S31_BURN, "sigma_isp": v["sigma_isp_s"],
            "sigma_b1": v["sigma_b1_s"], "corr": v["corr"], "k_wind": w["k_wind"],
            "sigma_k_wind": (k_at(55.0) - k_at(45.0)) / (2 * 1.96)}


def model_form_sigma_upper(mf_sigma):
    """Model discrepancy at the one-sided upper 95 % bound of the held-out apogee model error
    sigma_d (profile likelihood with r_i ~ N(0, sigma_i^2 + sigma_d^2), without the common
    calibration term; results/vsb30_zstats.json), as a common impulse factor scaled from the
    apogee-equivalent factor of model_form_sigma()."""
    v = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    z = json.load(open(os.path.join(ROOT, "results", "vsb30_zstats.json")))
    rms = v["stats_held_out"]["apogee_model_pct"]["rms"]
    up = z["sigma_d_upper95_pct"]
    return {"rms_pct": rms, "upper_pct": up, "upper_with_cal_pct": z["sigma_d_upper95_with_cal_pct"],
            "flight_sigma_pct": z["mean_flight_sigma_pct"], "sigma_impulse": mf_sigma * up / rms}


def cal_loo():
    """Booster parameters sampled with the leave-one-out spread of the single-flight calibration
    (run_vsb30_validation.py), centred on the V01 calibration; k_w as in calibration_uncertainty()."""
    v = json.load(open(os.path.join(ROOT, "results", "vsb30_validation.json")))
    rows = [x for x in v["loo"] if "isp" in x and "b1" in x]
    import statistics as st
    i_, b_ = [x["isp"] for x in rows], [x["b1"] for x in rows]
    mi, mb = st.mean(i_), st.mean(b_)
    si, sb = st.stdev(i_), st.stdev(b_)
    corr = sum((a - mi) * (b - mb) for a, b in zip(i_, b_)) / ((len(i_) - 1) * si * sb)
    cu = calibration_uncertainty()
    return dict(cu, sigma_isp=si, sigma_b1=sb, corr=corr)


def cal_multiflight():
    """Booster parameters sampled from the multi-flight calibration (run_multiflight_calibration.py)."""
    m = json.load(open(os.path.join(ROOT, "results", "multiflight_calibration.json")))
    cu = calibration_uncertainty()
    return dict(cu, isp=m["isp"], b1=m["b1"], sigma_isp=m["sigma_isp"], sigma_b1=m["sigma_b1"], corr=m["corr"])


def paired(base, other):
    """Paired comparison of the first min(n) runs (common random numbers): difference of the
    probabilities and exact McNemar p-value, for the window and for tau_use >= 3 s with the link."""
    from scipy.stats import binomtest
    out = {}
    for key in ("met", "use_ok"):
        a, b = base[key], other[key]
        n = min(len(a), len(b))
        n10 = sum(1 for i in range(n) if a[i] == "1" and b[i] == "0")
        n01 = sum(1 for i in range(n) if a[i] == "0" and b[i] == "1")
        p = binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 else 1.0
        out[key] = {"n": n, "diff": (n01 - n10) / n, "n10": n10, "n01": n01, "p_mcnemar": float(p)}
    return out


def run(concepts, n, k_wind, scale=1.0, seed=SEED, tag="", mf_sigma=0.0, cal=None, beta_scale=1.0,
        msis=False, att_deg=0.0):
    res = {}
    for c in concepts:
        # common random numbers: every concept and every set restarts the same
        # streams, so run i of any set uses the aleatory draws of run i of the base set
        rng = random.Random(seed)
        rng_mf = random.Random(seed + 1000)
        rng_cal = random.Random(seed + 2000)
        rng_att = random.Random(seed + 3000)
        draws = [D.sample(rng, k_wind, scale) for _ in range(n)]
        if beta_scale != 1.0:   # wider test-vehicle ballistic-coefficient dispersion
            for d in draws:
                d["beta"] = max(1.0 + (d["beta"] - 1.0) * beta_scale, 0.3)
        if cal is not None:     # epistemic: calibration and k_w uncertainty, one world per run
            for d in draws:
                z1, z2 = rng_cal.gauss(0, 1), rng_cal.gauss(0, 1)
                d["isp"] = cal["isp"] + cal["sigma_isp"] * z1
                # equivalent S31 burn time kept within the bounds of the calibration (4-20 s), as in shyma.calibration
                d["b1"] = min(max(cal["b1"] + cal["sigma_b1"] * (cal["corr"] * z1 + math.sqrt(1 - cal["corr"] ** 2) * z2),
                                  4.0), 20.0)
                d["k_wind"] = max(cal["k_wind"] + cal["sigma_k_wind"] * rng_cal.gauss(0, 1), 0.0)
        if mf_sigma > 0.0:      # model-form error: common impulse factor, independent stream
            for d in draws:
                f = rng_mf.gauss(1.0, mf_sigma)
                d["imp1"] *= f
                d["imp2"] *= f
        if msis:
            for d in draws:
                d["msis"] = True
        if att_deg > 0.0:       # attitude error at S30 ignition: two Gaussian components per axis
            for d in draws:
                a1, a2 = rng_att.gauss(0, att_deg), rng_att.gauss(0, att_deg)
                d["att"] = (math.hypot(a1, a2), math.atan2(a2, a1))
        with Pool(int(os.environ.get("MC_PROCS", "2"))) as p:
            out = p.map(one, [(c, d) for d in draws])
        res[c] = analyse(out)
        print(tag, c, {k: v for k, v in res[c].items() if k not in ("xy", "window_samples", "failures")},
              flush=True)
    return res


SETS = ("with_wind", "no_wind", "doubled_sigma", "model_form", "calibration", "beta_25",
        "msis", "model_form_upper", "calibration_loo", "calibration_multiflight", "ign_att_1deg", "ign_att_2deg")

if __name__ == "__main__":
    t0 = time.time()
    k = k_wind_calibrated()
    # MC_SETS=a,b,... recomputes only those sets and keeps the others of results/montecarlo.json
    only = [x for x in os.environ.get("MC_SETS", "").split(",") if x]
    path = os.path.join(ROOT, "results", "montecarlo.json")
    res = json.load(open(path)) if only else {}
    res.update({"N": N, "sigma": D.SIG, "k_wind": k, "seed": SEED, "ant_loss_db": ANT_LOSS_DB,
                "fade_db": FADE_DB, "nominal_elevation": cs.NOMINAL})
    todo = only or SETS
    if "with_wind" in todo:
        res["with_wind"] = run("ABCD", N, k, tag="wind")
    if "no_wind" in todo:
        res["no_wind"] = run("ABCD", N // 2, 0.0, seed=SEED, tag="nowind")
    if "doubled_sigma" in todo:
        res["doubled_sigma"] = run("ABCD", N_2X, k, scale=2.0, seed=SEED, tag="2x")
    res["model_form_sigma_impulse"] = model_form_sigma()
    if "model_form" in todo:
        res["model_form"] = run("ABCD", N // 2, k, seed=SEED, tag="mf", mf_sigma=res["model_form_sigma_impulse"])
    res["calibration_uncertainty"] = calibration_uncertainty()
    if "calibration" in todo:
        res["calibration"] = run("ABCD", N // 2, k, seed=SEED, tag="cal", cal=res["calibration_uncertainty"])
    if "beta_25" in todo:
        res["beta_25"] = run("AD", N // 2, k, seed=SEED, tag="beta25", beta_scale=2.5)
    if "msis" in todo:
        res["msis"] = run("ABCD", N // 2, k, seed=SEED, tag="msis", msis=True)
    res["model_form_upper_sigma"] = model_form_sigma_upper(res["model_form_sigma_impulse"])
    if "model_form_upper" in todo:
        res["model_form_upper"] = run("ABCD", N // 2, k, seed=SEED, tag="mfu",
                                      mf_sigma=res["model_form_upper_sigma"]["sigma_impulse"])
    res["calibration_loo_uncertainty"] = cal_loo()
    if "calibration_loo" in todo:
        res["calibration_loo"] = run("ABCD", N // 2, k, seed=SEED, tag="calloo", cal=res["calibration_loo_uncertainty"])
    res["calibration_multiflight_uncertainty"] = cal_multiflight()
    if "calibration_multiflight" in todo:
        # fitted with the NRLMSISE-00 site atmospheres, so flown with the NRLMSISE-00 atmosphere of CLA
        res["calibration_multiflight"] = run("ABCD", N // 2, k, seed=SEED, tag="calmf",
                                             cal=res["calibration_multiflight_uncertainty"], msis=True)
    # nominal state of the stack at S30 ignition (Concept D)
    v_, tr_ = cs.fly("D", 400.0, cs.NOMINAL["D"], rec=1)
    ti = v_.stages[1].t_ign
    ki = min(range(len(tr_.t)), key=lambda j: abs(tr_.t[j] - ti))
    res["d_ignition_state"] = {"t_s": ti, "alt_km": tr_.alt[ki] / 1e3, "mach": tr_.mach[ki],
                               "q_kpa": tr_.q[ki] / 1e3, "alpha_deg": tr_.alpha[ki]}
    if "ign_att_1deg" in todo:
        res["ign_att_1deg"] = run("D", N // 2, k, seed=SEED, tag="att1", att_deg=1.0)
    if "ign_att_2deg" in todo:
        res["ign_att_2deg"] = run("D", N // 2, k, seed=SEED, tag="att2", att_deg=2.0)
    res["paired_vs_base"] = {key: {c: paired(res["with_wind"][c], res[key][c]) for c in res[key]}
                             for key in SETS[1:] if key in res and key != "no_wind"}
    for key in SETS[1:]:
        for c in res.get(key, {}):
            res[key][c].pop("xy", None)
            res[key][c].pop("window_samples", None)
    res["cpu_s" if not only else "cpu_s_partial"] = time.time() - t0
    with open(os.path.join(ROOT, "results", "montecarlo.json"), "w") as f:
        json.dump(res, f)
    print("wall", res["cpu_s"])
