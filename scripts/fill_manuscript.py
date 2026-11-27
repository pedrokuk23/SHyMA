"""Fill the @@KEY@@ placeholders of paper/main.tpl.tex and paper/supplement.tpl.tex from results/*.json and
write paper/main.tex (and the same for paper_acta/, the Acta Astronautica version). Every number of the manuscript that comes from a
computation is produced here (results/manuscript_values.json lists them). In the public repository the
manuscript templates are absent: the script then writes only results/manuscript_values.json, which maps every
placeholder name of the article to its value.
"""
import json
import math
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)


def load(name):
    with open(os.path.join(ROOT, "results", name + ".json")) as f:
        return json.load(f)


def f(x, n=1):
    return ("%%.%df" % n) % x


def sgn(x, n=1):
    s = f(x, n)
    return ("$+$" + s) if x > 0 else ("$-$" + s.lstrip("-") if x < 0 else s)


def pct_ci(lo, hi, n=3):
    return "[%s, %s]" % (f(lo, n), f(hi, n))


def _new_round4_keys(K, aero):
    """Aerodynamic representation, mass margins and selection under uncertainty (round 4)."""
    # aerodynamic representation: frozen vs recalibrated propulsion (run_aero_sensitivity.py)
    nm = {"cn_alpha": "$C_{N\\alpha}$", "margin_cal": "$\\ell$", "cmq": "$C_{mq}$", "cd": "$C_D(M)$"}
    rows = []
    for k in ("nominal", "cn_alpha_x0.5", "cn_alpha_x2.0", "margin_cal_x0.5", "margin_cal_x2.0", "cmq_x0.5",
              "cmq_x2.0", "cd_x0.8", "cd_x1.2"):
        r = aero[k]
        lab = "Generic (nominal)" if k == "nominal" else "%s $\\times$%s" % (nm[k.rsplit("_x", 1)[0]], k.rsplit("_x", 1)[1])
        fr = r.get("frozen", r)
        rows.append("%s & %s & %s & %s & %s & %s & %s \\\\" % (
            lab, sgn(r["held_out"]["mean_pct"], 1), sgn(fr["held_out"]["mean_pct"], 1),
            f(r["dwell_A"], 2), f(fr["dwell_A"], 2), f(r["dwell_D"], 2), f(fr["dwell_D"], 2)))
    ga = load("geometric_aero")
    rows.append("\\midrule")
    for k, lab in (("geometric, 3 fins", "Geometric, 3 fins"), ("geometric, 4 fins", "Geometric, 4 fins"),
                   ("geometric, fins x0.9", "Geometric, fins $\\times$0.9"),
                   ("geometric, fins x1.1", "Geometric, fins $\\times$1.1")):
        r = ga["variants"][k]
        rows.append("%s & %s & -- & %s & -- & %s & -- \\\\" % (lab, sgn(r["held_out"]["mean_pct"], 1),
                                                             f(r["dwell_A"], 2), f(r["dwell_D"], 2)))
    gd = load("geometric_drag")
    for k, lab in (("geometric drag", "Geometric drag"), ("geometric drag, t/c 0.05", "Geometric drag, $t/c = 0.05$"),
                   ("geometric drag, exit ratio 0.5", "Geometric drag, $A_e/A_b = 0.5$"),
                   ("geometric drag and stability", "Geometric drag and derivatives")):
        r = gd["variants"][k]
        rows.append("%s & %s & -- & %s & -- & %s & -- \\\\" % (lab, sgn(r["held_out"]["mean_pct"], 1),
                                                             f(r["dwell_A"], 2), f(r["dwell_D"], 2)))
    K["AE_ROWS"] = "\n".join(rows)
    best = gd["variants"]["geometric drag and stability"]
    K["GD_ISP"], K["GD_ISP0"] = f(best["isp"], 1), f(load("aero_sensitivity")["nominal"]["isp"], 1)
    K["GD_HO_RMS"] = f(best["held_out"]["rms_pct"], 1)
    K["GD_HO_MEAN"] = sgn(best["held_out"]["mean_pct"], 1)
    K["GD_A"], K["GD_D"] = f(best["dwell_A"], 1), f(best["dwell_D"], 1)
    K["GD_AR"], K["GD_DR"] = f(best["A_retuned"]["dwell_s"], 1), f(best["D_retuned"]["dwell_s"], 1)
    K["GD_A_LO"], K["GD_A_HI"] = f(gd["dwell_range_A"][0], 1), f(gd["dwell_range_A"][1], 1)
    K["GD_D_LO"], K["GD_D_HI"] = f(gd["dwell_range_D"][0], 1), f(gd["dwell_range_D"][1], 1)
    K["GD_AR_LO"], K["GD_AR_HI"] = f(gd["dwell_retuned_range_A"][0], 1), f(gd["dwell_retuned_range_A"][1], 1)
    K["GD_DR_LO"], K["GD_DR_HI"] = f(gd["dwell_retuned_range_D"][0], 1), f(gd["dwell_retuned_range_D"][1], 1)
    K["GD_EL_A"] = "%.0f" % best["A_retuned"]["el"]
    K["GD_EL_D"], K["GD_CO_D"] = "%.0f" % best["D_retuned"]["el"], "%.0f" % best["D_retuned"]["coast"]
    tg = {round(x["mach"], 2): x for x in gd["table"]}
    K["GD_CD03"], K["GD_CD03ON"] = f(tg[0.3]["full_off"]["total"], 2), f(tg[0.3]["full_on"], 2)
    K["GD_CD2"], K["GD_CD6U"] = f(tg[2.0]["full_off"]["total"], 2), f(tg[6.0]["upper_off"], 2)
    tl, tb = ga["tables"]["full, lift-off"]["rows"], ga["tables"]["upper, S30 burnout"]["rows"]
    tu = ga["tables"]["upper, S30 ignition"]["rows"]
    by = lambda rows, m: [x for x in rows if abs(x["mach"] - m) < 1e-9][0]  # noqa: E731
    K["GA_CNA_03"], K["GA_CNA_7"] = f(by(tl, 0.3)["cn_alpha"], 0), f(by(tl, 7.0)["cn_alpha"], 0)
    K["GA_CNAU_03"], K["GA_CNAU_7"] = f(by(tu, 0.3)["cn_alpha"], 0), f(by(tu, 7.0)["cn_alpha"], 1)
    K["GA_SM_L03"], K["GA_SM_L3"] = f(by(tl, 0.3)["margin_cal"], 1), f(by(tl, 3.0)["margin_cal"], 1)
    K["GA_SM_UI6"], K["GA_SM_UB6"] = f(by(tu, 6.0)["margin_cal"], 1), f(by(tb, 6.0)["margin_cal"], 1)
    K["GA_CMQ_LO"] = sgn(max(x["cmq"] for t_ in (tl, tb) for x in t_), 0)
    K["GA_CMQ_HI"] = sgn(min(x["cmq"] for t_ in (tl, tb) for x in t_), 0)
    K["GA_DA_LO"], K["GA_DA_HI"] = f(ga["dwell_range_A"][0], 2), f(ga["dwell_range_A"][1], 2)
    K["GA_DD_LO"], K["GA_DD_HI"] = f(ga["dwell_range_D"][0], 2), f(ga["dwell_range_D"][1], 2)
    gu = [ga["variants"][k]["uwe_km_per_mps"] / ga["variants"]["generic"]["uwe_km_per_mps"] - 1
          for k in ga["variants"] if k != "generic"]
    K["GA_UWE_LO"], K["GA_UWE_HI"] = "%.0f" % (100 * min(gu)), "%.0f" % (100 * max(gu))
    sc = load("stability_check")["inverse"]
    K["GA_INV_ERR"] = "%.2f" % max(abs(sc["err_cm_alpha_pct"]), abs(sc["err_cmq_pct"]))
    for cc in "AD":
        for gk, m in (("STR", "stab_recal"), ("STF", "stab_frozen"), ("CDR", "cd_recal"), ("CDF", "cd_frozen"),
                     ("CDT", "cd_retuned")):
            lo_, hi_ = aero["dwell_%s_%s" % (m, cc)] if m != "cd_retuned" else aero["dwell_cd_retuned_" + cc]
            K["AE_%s_%s_LO" % (gk, cc)], K["AE_%s_%s_HI" % (gk, cc)] = f(lo_, 1), f(hi_, 1)
    K["AE_HO_NOM"] = sgn(aero["nominal"]["held_out"]["mean_pct"], 1)
    hr = [aero[k]["held_out"]["mean_pct"] for k in ("cd_x0.8", "cd_x1.2")]
    hf = [aero[k]["frozen"]["held_out"]["mean_pct"] for k in ("cd_x0.8", "cd_x1.2")]
    K["AE_HO_CDR_LO"], K["AE_HO_CDR_HI"] = sgn(min(hr), 1), sgn(max(hr), 1)
    K["AE_HO_CDF_LO"], K["AE_HO_CDF_HI"] = sgn(min(hf), 1), sgn(max(hf), 1)
    # mass margins of the architectures (run_mass_budget.py)
    mb = load("mass_budget")["concepts"]
    for cc in "ABCD":
        fb = mb[cc]["fixed-beta"]
        K["MB_CW_" + cc] = "%.0f" % math.floor(fb["ceiling_window_kg"])
        K["MB_CU_" + cc] = "%.0f" % math.floor(fb["ceiling_useful_kg"])
        K["MB_FU_" + cc] = "%.0f" % math.ceil(fb["floor_useful_kg"])
    for cc in "AB":
        K["MB_CG_" + cc] = "%.0f" % math.floor(mb[cc]["fixed-geometry"]["ceiling_window_kg"])
    # selection under uncertainty (run_robust_selection.py)
    rs = load("robust_selection")
    K["RS_NS"], K["RS_NC"] = "%d" % rs["n_screen"], "%d" % rs["n_confirm"]
    K["RS_NCAND"] = "%d" % sum(len(r["screen"]) for r in rs["concepts"].values())
    K["RS_WALL"] = "%.0f" % (rs["cpu_s"] / 60)
    K["RS_PMIN"] = f(rs["p_min"], 2)

    def xs(c, x):
        return ("%.0f\\textdegree" % x["el"]) + (", %.0f~s" % x["coast"] if c == "D" else "")
    rows = []
    for c in "ABCD":
        r = rs["concepts"][c]
        seen = []
        for tag, lab in (("nominal", "nominal"), ("primary", "primary"), ("R1", "max.\\ $P$"), ("R2", "min.\\ range"), ("R3", "max.\\ 5th pct.")):
            if tag not in r["confirm"]:
                continue
            v = r["confirm"][tag]
            if v["x"] in seen:
                continue
            seen.append(v["x"])
            labs = [l for t, l in (("nominal", "nominal"), ("primary", "primary"), ("R1", "max.\\ $P$"), ("R2", "min.\\ range"),
                                   ("R3", "max.\\ 5th pct.")) if t in r["confirm"] and r["confirm"][t]["x"] == v["x"]]
            rows.append("%s & %s & %s & %s [%s, %s] & %s & %s & %s & %.0f \\\\" % (
                c, ", ".join(labs), xs(c, v["x"]), f(v["p_req"], 3), f(v["p_req_ci95"][0], 3), f(v["p_req_ci95"][1], 3),
                f(v["p_req_rec"], 3), f(v["tau_use_p05"], 1), f(v["tau_use_p50"], 1), v["range_mean_km"]))
    K["RS_ROWS"] = "\n".join(rows)
    for c in "ABCD":
        cf = rs["concepts"][c]["confirm"]
        K["RS_P_" + c] = f(cf["nominal"]["p_req"], 3)
        K["RS_P05_" + c] = f(cf["nominal"]["tau_use_p05"], 1)
        K["RS_P50_" + c] = f(cf["nominal"]["tau_use_p50"], 1)
        K["RS_R_" + c] = "%.0f" % cf["nominal"]["range_mean_km"]
    a3 = rs["concepts"]["A"]["confirm"]["R3"]
    K["RS_A3_EL"], K["RS_A3_P"], K["RS_A3_P05"] = "%.0f" % a3["x"]["el"], f(a3["p_req"], 3), f(a3["tau_use_p05"], 1)
    K["RS_P05_A2"] = f(rs["concepts"]["A"]["confirm"]["nominal"]["tau_use_p05"], 2)
    d3 = rs["concepts"]["D"]["confirm"]["R3"]
    nf = round((1.0 - d3["p_req"]) * d3["n"])
    K["RS_D3_FAIL"] = ("no failure" if nf == 0 else ("%d failure%s" % (nf, "" if nf == 1 else "s")))
    K["RS_D3_EL"], K["RS_D3_CO"] = "%.0f" % d3["x"]["el"], "%.0f" % d3["x"]["coast"]
    K["RS_D3_P"], K["RS_D3_PLO"] = f(d3["p_req"], 3), f(d3["p_req_ci95"][0], 3)
    K["RS_D3_P05"], K["RS_D3_P50"], K["RS_D3_R"] = f(d3["tau_use_p05"], 1), f(d3["tau_use_p50"], 1), "%.0f" % d3["range_mean_km"]
    for c in "BC":
        sc = [q for q in rs["concepts"][c]["screen"] if q["x"]["el"] < 84.0]
        K["RS_LINK0_" + c] = "%.0f--%.0f" % (min(q["x"]["el"] for q in sc), max(q["x"]["el"] for q in sc))
        assert max(q["p_req"] for q in sc) == 0.0
        rec_ok = [q for q in sc if q["p_req_rec"] >= 0.95]
        K["RS_REC_" + c] = "%.0f" % min(q["x"]["el"] for q in rec_ok)
    if K.get("RS_REC_B") == K.get("RS_REC_C"):
        K["RS_REC_TXT"] = "from an elevation of %s\\textdegree\\ for both" % K["RS_REC_B"]
    else:
        K["RS_REC_TXT"] = "from an elevation of %s\\textdegree\\ for B and of %s\\textdegree\\ for C" % (K.get("RS_REC_B"), K.get("RS_REC_C"))


def main():
    K = {}
    ver, nesc, jsb = load("verification"), load("nesc_checkcases"), load("jsbsim_verification")
    vsb, wnd, hy = load("vsb30_validation"), load("wind_calibration"), load("hyshot_validation")
    cs, aero = load("casestudy"), load("aero_sensitivity")
    cal = json.load(open(os.path.join(ROOT, "shyma", "calibration.json")))

    # calibration
    K["ISP"] = f(cal["isp_vsb30"], 1)
    _d = [100 * (1 - cal["isp_vsb30"] / x) for x in (260.0, 263.0)]   # IAE vacuum Isp of S31 and S30 (Palmerio, 2017)
    K["ISP_IAE_D"] = ("%.1f--%.1f\\%% below" % (min(_d), max(_d))) if min(_d) > 0 else ("within %.1f\\%% of" % max(abs(x) for x in _d))
    K["TB1"] = f(cal["s31_burn"], 2)
    J = cal["jacobian_km_per_s"]
    K["JAC"] = r"$\begin{bmatrix} %s & %s \\ %s & %s \end{bmatrix}$" % (
        f(J[0][0], 2), sgn(J[0][1], 2).replace("$", ""), f(J[1][0], 2), f(J[1][1], 2))
    cu = vsb["calibration_uncertainty"]
    K["SIG_ISP"] = f(cu["sigma_isp_s"], 1)
    K["SIG_TB1"] = f(cu["sigma_b1_s"], 1)
    K["SIG_TB1_NP"] = f(cu["without_residual_profile"]["sigma_b1_s"], 1) if "without_residual_profile" in cu else "--"
    K["SIG_ISP_NP"] = f(cu["without_residual_profile"]["sigma_isp_s"], 1) if "without_residual_profile" in cu else "--"
    K["CAL_RHO"] = f(cu.get("obs_corr", 0.0), 2)

    # verification
    def e(x):
        m, ex = ("%.1e" % x).split("e")
        return r"$%s\times10^{%d}$" % (m, int(ex))
    K["KEP_E20"] = e(ver["kepler_dt_20.00"]["pos_err_m_t200"])
    K["KEP_E25"] = e(ver["kepler_dt_2.50"]["pos_err_m_t200"])
    K["KEP_ORD"] = f(sum(ver["kepler_observed_order"]) / len(ver["kepler_observed_order"]), 2)
    K["KEP_ORD_LO"], K["KEP_ORD_HI"] = f(min(ver["kepler_observed_order"]), 2), f(max(ver["kepler_observed_order"]), 2)
    imp = max(ver[k]["impact_diff_m"] for k in ver if k.startswith(("vsb30_dt", "fly5_dt")))
    apo = max(abs(ver[k]["apogee_diff_m"]) for k in ver if k.startswith(("vsb30_dt", "fly5_dt")))
    K["VER_IMP"] = f(math.ceil(imp * 10) / 10, 1)
    K["VER_APO"] = f(math.ceil(apo * 100) / 100, 2)
    win = max(abs(ver["window_%s" % c]["0.050"] - ver["window_%s" % c]["0.005"]) for c in "ABCD" if "window_%s" % c in ver)
    K["VER_WIN"] = "%d" % max(1, math.ceil(win * 1000))

    # NESC
    K["N09POS"] = f(nesc["case09"]["pos_err_final_m"], 2)
    K["N10POS"] = f(nesc["case10"]["pos_err_final_m"], 2)
    m17 = nesc["case17"]["metrics"]
    alt17 = nesc["case17"]["checkpoints"]["134.17"]["alt_m"]
    K["N17ALT"] = "%d" % round(abs(m17["alt_m"]["final_err"]))
    K["N17ALTPCT"] = f(100 * abs(m17["alt_m"]["final_err"]) / alt17["ref_median"], 2)
    K["N17ALTSP"] = "%d" % round(m17["alt_m"]["max_spread"])
    K["N17VEL"] = f(max(m17[q]["max_abs_err"] for q in ("vn_mps", "ve_mps", "vd_mps")), 1)
    K["N17VELSP"] = f(max(m17[q]["max_spread"] for q in ("vn_mps", "ve_mps", "vd_mps")), 1)
    K["N17ALTX"] = f(m17["alt_m"]["max_abs_err"] / m17["alt_m"]["max_spread"], 1)
    K["N17VELX"] = f(max(m17[q]["max_abs_err"] / m17[q]["max_spread"] for q in ("ve_mps", "vd_mps")), 1)
    K["N17LON"] = "%.4f" % m17["lon_deg"]["max_abs_err"]
    K["N17LONX"] = f(m17["lon_deg"]["max_abs_err"] / m17["lon_deg"]["max_spread"], 0)
    K["N17LONSP"] = "%.5f" % m17["lon_deg"]["max_spread"]

    def unit_fmt(q, v):
        if q == "alt_m":
            return ("%.3f m" % v) if v < 0.1 else ("%.2f m" % v if v < 10 else "%.0f m" % v)
        if q.endswith("_mps"):
            return ("%.4f m/s" % v) if v < 0.1 else "%.2f m/s" % v
        if v < 1e-3:
            m, ex = ("%.1e" % v).split("e")
            return r"$%s\times10^{%d}$" % (m, int(ex))
        return "%.4f" % v
    names = {"alt_m": "altitude", "ve_mps": "$v_E$", "vd_mps": "$v_D$", "vn_mps": "$v_N$", "mach": "Mach",
             "lon_deg": "longitude (deg)"}
    rows = []
    for key, lab, qs in (("case05", "5", ("alt_m", "vd_mps", "mach")),
                         ("case09", "9", ("alt_m", "ve_mps", "vd_mps", "mach")),
                         ("case10", "10", ("alt_m", "vn_mps", "vd_mps")),
                         ("case17", "17", ("alt_m", "ve_mps", "vd_mps", "lon_deg", "mach"))):
        c = nesc[key]
        for q in qs:
            mm = c["metrics"][q]
            rows.append("Case %s, %s & %s & %s & %d \\\\" % (lab, names[q], unit_fmt(q, mm["max_abs_err"]),
                                                           unit_fmt(q, mm["max_spread"]), mm["n_ref"]))
        if key != "case17":
            rows.append("\\midrule")
    K["NESC_ROWS"] = "\n".join(rows)

    # JSBSim
    K["JSB_W1"] = f(jsb["entry85_wgs84"]["impact_sep_km"] * 1e3, 1)
    K["JSB_W2"] = f(jsb["apogee250_wgs84"]["impact_sep_km"] * 1e3, 1)
    K["JSB_S1"] = f(jsb["entry85"]["impact_sep_km"], 2)
    K["JSB_S2"] = f(jsb["apogee250"]["impact_sep_km"], 2)

    # VSB-30
    st = vsb["stats_held_out"]
    a, r, t = st["apogee_model_pct"], st["range_model_km"], st["t100_model_pct"]
    K.update({"V_AM": sgn(a["mean"]), "V_AR": f(a["rms"]), "V_AX": f(a["max_abs"]),
              "V_RM": sgn(r["mean"]), "V_RR": f(r["rms"]), "V_RX": f(r["max_abs"]),
              "V_TM": sgn(t["mean"]), "V_TR": f(t["rms"]), "V_TX": f(t["max_abs"])})
    K["V_AR"] = f(a["rms"])
    ts = vsb["tests"]
    K["V_DCI0"] = sgn(ts["rms_diff_ci95"][0])
    K["V_DCI1"] = sgn(ts["rms_diff_ci95"][1])
    K["V_WP"] = f(ts["wilcoxon_abs_p"], 2)
    K["V_CLA"] = sgn(st["apogee_model_pct_CLA"]["mean"])
    K["V_ESR"] = sgn(st["apogee_model_pct_Esrange"]["mean"])
    K["V_SPH"] = f(vsb["sphere"]["stats_held_out"]["apogee_model_pct"]["rms"])
    loo = vsb["loo"]
    ar = sorted(x["apogee_rms_pct"] for x in loo)
    K["LOO_AMIN"], K["LOO_AMAX"] = f(ar[0]), f(ar[-1])
    K["LOO_AMED"] = f(0.5 * (ar[4] + ar[5]))
    K["LOO_IMIN"] = "%.0f" % math.floor(min(x["isp"] for x in loo))
    K["LOO_IMAX"] = "%.0f" % math.ceil(max(x["isp"] for x in loo))
    K["LOO_BMIN"] = f(min(x["b1"] for x in loo))
    K["LOO_BMAX"] = f(max(x["b1"] for x in loo))
    K["LOO_RMIN"] = "%.0f" % min(x["range_rms_km"] for x in loo)
    K["LOO_RMAX"] = "%.0f" % max(x["range_rms_km"] for x in loo)
    nz = vsb["normalised_residuals"]
    K["Z_A"], K["Z_R"] = f(nz["rms_z_apogee"], 2), f(nz["rms_z_range"], 2)
    K["Z_PA"], K["Z_PR"] = f(nz["p_apogee"], 3), f(nz["p_range"], 2)
    K["Z_PR2"] = f(2 * min(nz["p_range"], 1 - nz["p_range"]), 2)   # two-sided chi^2 p of the range
    zs = load("vsb30_zstats")
    K["ZS_CAL"] = f(zs["sigma_cal_apogee_km"])
    K["ZS_P"] = f(zs["p"], 2)

    # wind
    K["KW"] = f(wnd["k_wind"], 2)
    K["KW_LO"], K["KW_HI"] = f(wnd["k_wind_ci95"][0], 2), f(wnd["k_wind_ci95"][1], 2)
    tab = {float(k): v for k, v in wnd["table"].items()}
    ks = sorted(tab)

    def k_at(target):
        for a0, b0 in zip(ks[:-1], ks[1:]):
            ra, rb = tab[a0]["r3_gauss"], tab[b0]["r3_gauss"]
            if (ra - target) * (rb - target) <= 0:
                return a0 + (target - ra) * (b0 - a0) / (rb - ra)
        return float("nan")
    K["KW_T45"], K["KW_T55"] = f(k_at(45.0), 2), f(k_at(55.0), 2)
    K["KW_EMP"] = f(wnd["k_wind_empirical_estimator"], 2)
    K["R3NW"] = f(wnd["no_wind"]["r3_gauss"])
    K["RSD"] = f(wnd["at_k"]["range_sd"])
    K["UWE"] = f(aero["nominal"]["uwe_km_per_mps"])
    K["UWE_LO"] = sgn(aero["uwe_change_pct_range"][0], 0)
    K["UWE_HI"] = sgn(aero["uwe_change_pct_range"][1], 0)
    d = wnd["weighting_demo"]
    K["WW_SHIFT"] = "%.0f" % d["unweighted_shift_km"]
    K["WW_IT"] = "%d" % d["iterations"]
    K["WW_EL"] = f(d["weighted_elevation"], 1)
    K["WW_AZ"] = f(d["weighted_azimuth"], 1)
    K["WW_MISS"] = "%.0f" % max(1, math.ceil(d["residual_miss_m"]))

    # HyShot
    n = hy["nominal"]
    K["HY_M35"], K["HY_M29"] = f(n["mach_35"], 2), f(n["mach_29"], 2)
    K["HY_EM35"] = f(100 * (n["mach_35"] - 7.75) / 7.75)
    K["HY_EM29"] = f(100 * (n["mach_29"] - 8.0) / 8.0)
    K["HY_DM35"], K["HY_DM29"] = f(n["mach_35"] - 7.75, 2), f(n["mach_29"] - 8.0, 2)
    K["HY_Q35"], K["HY_Q29"] = f(n["q_35_kpa"]), f(n["q_29_kpa"])
    K["HY_EQ35"] = "%.0f" % (100 * (n["q_35_kpa"] - 21.0) / 21.0)
    K["HY_EQ29"] = "%.0f" % (100 * (n["q_29_kpa"] - 58.5) / 58.5)
    K["HY_TOPEN"] = f(n["t_35km"])
    dto = n["t_35km"] - 537.5
    K["HY_DTO"] = f(abs(dto))
    K["HY_DTOS"] = sgn(dto)
    K["HY_DTOWORD"] = "late" if dto >= 0 else "early"
    K["HY_DTO_SE"] = f(hy["sensitivity"]["space_end_state"]["t_35km"] - 537.5)
    K["HY_SE_DT"] = f(hy["space_end_state"]["t"] - hy["entry_state"]["t"])
    r5 = hy["ref_alt_2005"]
    K["HY5_DUR"] = f(r5["duration_35_29"], 2)
    K["HY5_EQ35"] = "%.0f" % (100 * (r5["q_35_kpa"] - 21.0) / 21.0)
    K["HY5_EQ29"] = "%.0f" % (100 * (r5["q_29_kpa"] - 58.5) / 58.5)
    K["HY5_EM35"] = f(100 * (r5["mach_35"] - 7.75) / 7.75)
    K["HY5_EM29"] = f(100 * (r5["mach_29"] - 8.0) / 8.0)
    K["HY_SE_M35"] = sgn(100 * (hy["sensitivity"]["space_end_state"]["mach_35"] - 7.75) / 7.75)
    K["HY_DUR"] = f(n["duration_35_29"])
    sb = [hy["sensitivity"][k] for k in ("beta_1500", "beta_5000", "beta_15000")]
    K["HY_MB0"], K["HY_MB1"] = f(min(x["mach_35"] for x in sb), 2), f(max(x["mach_35"] for x in sb), 2)
    K["HY_MC0"], K["HY_MC1"] = f(min(x["mach_29"] for x in sb), 2), f(max(x["mach_29"] for x in sb), 2)
    K["HY_QB0"], K["HY_QB1"] = f(min(x["q_35_kpa"] for x in sb)), f(max(x["q_35_kpa"] for x in sb))
    K["HY_QC0"], K["HY_QC1"] = f(min(x["q_29_kpa"] for x in sb)), f(max(x["q_29_kpa"] for x in sb))
    K["HY_MBLO"] = sgn(100 * (min(x["mach_35"] for x in sb) - 7.75) / 7.75)
    K["HY_MBHI"] = sgn(100 * (max(x["mach_35"] for x in sb) - 7.75) / 7.75)
    K["HY_MCLO"] = sgn(100 * (min(x["mach_29"] for x in sb) - 8.0) / 8.0)
    K["HY_DMMIN"] = sgn(min(x["mach_29"] for x in sb) - 8.0, 2)
    K["HY_MCHI"] = sgn(100 * (max(x["mach_29"] for x in sb) - 8.0) / 8.0)
    K["HY_EMAX"] = f(max([abs(100 * (x["mach_35"] - 7.75) / 7.75) for x in sb] +
                         [abs(100 * (x["mach_29"] - 8.0) / 8.0) for x in sb]), 1)

    # case study
    N = cs["nominal"]
    K["TBO"] = f(N["A"]["t_burnout"], 1)
    K["B84H"] = f(cs["A_at_84"]["burnout_alt_km"])
    K["B84M"] = f(cs["A_at_84"]["burnout_mach"], 2)
    K["A_BOH"] = f(N["A"]["burnout_alt_km"])
    K["A_BOM"] = f(N["A"]["burnout_mach"], 2)
    K["A_BOQ"] = "%.0f" % N["A"]["burnout_q_kpa"]
    K["A_SHR"] = f(N["A"]["qdot_shroud_nose_boost_MW"], 1)
    for P in "ABCD":
        x = N[P]
        K[P + "_TAU"] = f(x["window_s"], 2)
        K[P + "_OBJ"] = f(x["objective_s"], 1)
        K[P + "_STD"] = f(x["steady_s"], 1)
        K[P + "_H0"], K[P + "_H1"] = f(x["window_alt_km"][0]), f(x["window_alt_km"][1])
        K[P + "_M0"], K[P + "_M1"] = f(x["window_mach"][0], 2), f(x["window_mach"][1], 2)
        K[P + "_QM"] = "%.0f" % x["window_q_mean_kpa"]
        K[P + "_FPA"] = sgn(x["window_fpa_entry_deg"], 0)
        dr = x["window_downrange_km"]
        K[P + "_DR"] = "%.0f--%.0f" % (min(dr), max(dr))
        K[P + "_APO"] = f(x["apogee_km"])
        K[P + "_R"] = "%.0f" % x["impact_range_km"]
        K[P + "_QMAX"] = "%.0f" % x["q_max_ascent_kpa"]
        K[P + "_QA"] = "%.0f" % x["qalpha_max_kpa_deg"]
        ws = x["window_stations"]["CLA"]
        K[P + "_EL0"], K[P + "_EL1"] = f(ws["el_min"]), f(ws["el_max"])
        K[P + "_EL0R"], K[P + "_EL1R"] = "%.0f" % ws["el_min"], "%.0f" % ws["el_max"]
        K[P + "_LM"] = f(ws["margin_min_db"])
        K[P + "_LMR"] = "%.0f" % ws["margin_min_db"]
        K[P + "_LMNET"] = f(ws["margin_min_db"] - 16.0)
        K[P + "_QS"] = f(x["qdot_le_window_MW"], 1)
        K[P + "_QHW"] = f(x["qdot_le_hw_window_MW"], 2)
        K[P + "_TW"] = "%.0f" % x["T_le_window_K"]
        K[P + "_T0"] = "%.0f" % x["T0_window_K"]
        K[P + "_QL"] = f(x["Q_le_to_window_end_MJ"], 1)
        K[P + "_TAU_MS"] = f(cs["sensitivity"]["msis_" + P]["window_s"], 2)
        sg = x["stages"]
        K[P + "_S31"] = "%.0f" % sg["S31"]["range_km"]
        if "S30" in sg:
            K[P + "_S30"] = "%.0f" % sg["S30"]["range_km"]
    K["MMAX"] = f(max(max(N[P]["window_mach"]) for P in "ABCD"), 2)
    tws = [N[P]["T_le_window_K"] for P in "ABCD"]
    t0s = [N[P]["T0_window_K"] for P in "ABCD"]
    K["TW_MIN"], K["TW_MAX"] = "%.0f" % min(tws), "%.0f" % max(tws)
    K["T0_MIN"], K["T0_MAX"] = "%.0f" % min(t0s), "%.0f" % max(t0s)
    th = load("thermal_model_check")
    rt = [100 * (th[P]["T0_cpg"] / th[P]["T0_tp"] - 1) for P in "ABCD"]
    rq = [100 * (th[P]["qhw_cpg"] / th[P]["qhw_tp"] - 1) for P in "ABCD"]
    K["CPG_T0"] = "%.0f--%.0f" % (min(rt), max(rt))
    K["CPG_Q"] = "%.0f--%.0f" % (min(rq), max(rq))
    S = cs["sensitivity"]
    for dd in range(5):
        K["SD_%d" % dd] = f(S["dsep_%d" % dd]["window_s"], 2)
    for P in "ABCD":
        for b in (2500, 5000, 10000, 15000):
            K["SB_%s%d" % (P, b)] = f(S["beta_%s_%d" % (P, b)]["window_s"], 2)
    K["SB_ADROP"] = f(S["beta_A_5000"]["window_s"] - S["beta_A_2500"]["window_s"], 1)
    tb = S["toward_CLBI_B"]
    K["CLBI_AZ"] = "%.0f" % tb["azimuth"]
    K["CLBI_EL0"] = f(tb["stations"]["CLBI"]["el_min"])
    K["CLBI_EL1"] = f(tb["stations"]["CLBI"]["el_max"])
    K["CLBI_EL1N"] = f(-tb["stations"]["CLBI"]["el_max"])
    K["CLBI_TAU"] = f(tb["window_s"], 2)
    K["CLBI_R"] = "%.0f" % tb["impact_range_km"]
    K["CS_WALL"] = "%.0f" % cs["cpu_s"]["total"]
    K["CS_WALL_MIN"] = "%.0f" % (cs["cpu_s"]["total"] / 60.0)
    g = cs["grid"]
    feas = {}
    for P in "ABC":
        feas[P] = [(p, el, g[P][i][j]) for i, p in enumerate(g["payloads"])
                   for j, el in enumerate(g["elevations"]) if g[P][i][j]["window_met"]]
        K["G%s_N" % P] = "%d" % len(feas[P])
    elsA = sorted(set(el for _, el, _ in feas["A"]))
    K["GA_EL"] = "%d--%d\\textdegree" % (elsA[0], elsA[-1])
    K["GA_RMIN"] = "%.0f" % min(s["impact_range_km"] for _, _, s in feas["A"])
    pB = sorted(set(p for p, _, _ in feas["B"]))
    K["GB_PTOP"] = "%d" % pB[-1]
    K["GB_PMAX"] = "%d" % (pB[-2] if len(pB) > 1 else pB[-1])
    K["GB_ELTOP"] = "%d" % min(el for p, el, _ in feas["B"] if p == pB[-1])
    K["GMAX"] = "%d" % max(p for P in "ABC" for p, _, _ in feas[P])
    for P in "ABC":
        op = [p for i, p in enumerate(g["payloads"]) for j, el in enumerate(g["elevations"])
              if g[P][i][j]["objective_met"]]
        K["GO_" + P] = "%d" % max(op) if op else "---"
    # Concept D
    K["D_EL"] = "%g" % N["D"]["elevation"]
    K["A_EL"] = "%.0f" % N["A"]["elevation"]
    K["D_COAST"] = "%.0f" % cs["coast_D"]
    K["D_PWR"] = f(N["D"]["t_burnout"] - N["D"]["window_t"][0], 1)
    K["D_MPK"] = f(N["D"]["burnout_mach"], 2)          # peak window Mach of D (S30 burnout inside the window)
    for c in (2, 10, 16, 20, 24, 27, 30):
        K["SC_%d" % c] = f(cs["sensitivity"]["coast_D_%d" % c]["window_s"], 2)
    gd = g["D"]
    def dpts(p, key):
        rows = gd[str(p)]
        return [(rows[i][j][key], el, c, rows[i][j]) for i, c in enumerate(gd["coasts"])
                for j, el in enumerate(gd["elevations"])]
    p400 = dpts(400, "window_s")
    K["GD_N400"] = "%d" % sum(1 for x in p400 if x[3]["window_met"])
    best = max(p400)
    K["GD_MAX400"] = f(best[0], 1)
    K["GD_BEST400"] = "%d\\textdegree, %d~s" % (best[1], best[2])
    K["GD_O350"] = f(max(x[0] for x in dpts(350, "objective_s")), 1)
    o300 = [x for x in dpts(300, "objective_s") if x[3]["objective_met"]]
    K["GD_O300_N"] = "%d" % len(o300)
    b3 = max(o300)
    K["GD_O300"], K["GD_O300_EL"], K["GD_O300_C"] = f(b3[0], 1), "%d" % b3[1], "%d" % b3[2]
    K["GD_O300_R"] = "%.0f--%.0f" % (min(x[3]["impact_range_km"] for x in o300),
                                     max(x[3]["impact_range_km"] for x in o300))
    od = [p for p in gd["payloads"] if any(x[3]["objective_met"] for x in dpts(p, "objective_s"))]
    K["GO_D"] = "%d" % max(od) if od else "---"
    # inverse design (scripts/run_optimizer.py)
    op = load("optimizer")
    for c in "ABCD":
        m4, o3, mc_ = op["max_dwell_400"][c], op["objective_350"][c], op["mass_ceiling"][c]
        K["OP4_" + c] = f(m4.get("threshold_s_fine", m4["threshold_s"]), 2)
        K["OP4_EL_" + c] = f(m4["x"]["el"], 2)
        K["OP4_R_" + c] = "%.0f" % m4["impact_range_km"]
        K["OP3_" + c] = f(o3["objective_s"], 1) if o3["objective_met"] else "not met (%s~s)" % f(o3["objective_s"], 1)
        K["OP3_EL_" + c] = f(o3["x"]["el"], 2)
        K["OPM_" + c] = "%.0f" % mc_["ceiling_kg"]
        K["OPN_R_" + c] = "%.0f" % N[c]["impact_range_km"]
        K["OP3T_" + c] = f(o3["objective_s"], 1) if o3["objective_met"] else "--"
    K["OP4_C_D"] = f(op["max_dwell_400"]["D"]["x"]["coast"], 1)
    K["OP3_C_D"] = f(op["objective_350"]["D"]["x"]["coast"], 1)
    K["OP4_GAIN_D"] = "%.0f" % (100 * (N["D"]["window_s"] / best[0] - 1))   # nominal (= optimum) vs best grid point
    assert abs(N["D"]["elevation"] - op["max_dwell_400"]["D"]["x"]["el"]) < 1e-9 and \
        abs(cs["coast_D"] - op["max_dwell_400"]["D"]["x"]["coast"]) < 1e-9, "Concept D nominal is not the optimum"
    K["OP_NEVAL_D"] = "%d" % op["max_dwell_400"]["D"]["n_eval"]
    K["OP_WALL"] = "%.0f" % (op["cpu_s"] / 60)
    # validation domain (scripts/run_validation_domain.py)
    vd = load("validation_domain")
    lab = [("elevation", "Launch elevation (\\textdegree)", 1), ("payload_kg", "Payload or test vehicle (kg)", 0),
           ("q_max_ascent_kpa", "Max.\\ ascent $q$ (kPa)", 0), ("qalpha_max_kpa_deg", "Max.\\ $q\\alpha$, boost (kPa\\,deg)", 0),
           ("burnout_mach", "Burnout Mach number", 2), ("apogee_km", "Apogee (km)", 0)]
    rows = []
    for k, name, d in lab:
        lo_, hi_ = vd["flight_range"][k]
        cells = []
        for c in "ABCD":
            x = vd["concepts"][c][k]
            t = ("%." + str(d) + "f") % x
            cells.append(t if vd["inside"][c][k] else "\\textbf{%s}" % t)
        rows.append("%s & %s--%s & %s \\\\" % (name, ("%." + str(d) + "f") % lo_, ("%." + str(d) + "f") % hi_,
                                             " & ".join(cells)))
    K["VD_ROWS"] = "\n".join(rows)
    K["VD_EL_LO"], K["VD_EL_HI"] = f(vd["flight_range"]["elevation"][0], 1), f(vd["flight_range"]["elevation"][1], 1)
    K["VD_QA_HI"] = "%.0f" % vd["flight_range"]["qalpha_max_kpa_deg"][1]
    K["VD_Q_HI"] = "%.0f" % vd["flight_range"]["q_max_ascent_kpa"][1]
    K["VD_QA_A"] = "%.0f" % vd["concepts"]["A"]["qalpha_max_kpa_deg"]
    K["VD_QA_D"] = "%.0f" % vd["concepts"]["D"]["qalpha_max_kpa_deg"]
    K["VD_QAX_A"] = f(vd["concepts"]["A"]["qalpha_max_kpa_deg"] / vd["flight_range"]["qalpha_max_kpa_deg"][1], 1)
    # aerodynamic sensitivity of the dwell of the depressed concepts
    K["AS_A_LO"], K["AS_A_HI"] = f(aero["dwell_range_A"][0], 2), f(aero["dwell_range_A"][1], 2)
    K["AS_D_LO"], K["AS_D_HI"] = f(aero["dwell_range_D"][0], 2), f(aero["dwell_range_D"][1], 2)
    _new_round4_keys(K, aero)
    # Sutton-Graves vs Fay-Riddell
    at = load("aerothermal_check")
    K["AT_C_LO"], K["AT_C_HI"] = sgn(100 * (at["ratio_cold_range"][0] - 1), 1), sgn(100 * (at["ratio_cold_range"][1] - 1), 1)
    K["AT_H_LO"], K["AT_H_HI"] = sgn(100 * (at["ratio_hot_range"][0] - 1), 1), sgn(100 * (at["ratio_hot_range"][1] - 1), 1)
    K["AT_MAXABS"] = f(100 * max(abs(r - 1) for r in at["ratio_cold_range"] + at["ratio_hot_range"]), 1)
    # model-form Monte Carlo
    mcj = load("montecarlo")
    K["MF_SIG"] = f(100 * mcj["model_form_sigma_impulse"], 2)
    K["MF_N"] = "%d" % mcj["model_form"]["A"]["n"]
    for c in "ABCD":
        K["MF_P" + c] = f(mcj["model_form"][c]["p_window"], 3)
        K["MF_P05_" + c] = f(mcj["model_form"][c]["window_p05"], 2)
    for c in "AD":
        K["MF_D" + c] = f(mcj["with_wind"][c]["p_window"] - mcj["model_form"][c]["p_window"], 2)
    fb = [x["burnout_mach"] for x in mcj["with_wind"]["A"]["failures"] if x.get("burnout_mach")]
    K["MA_FAILM"] = f(round(max(fb), 1), 1) if fb else "---"
    K["M2D_OBJN"] = "%d" % round(mcj["doubled_sigma"]["D"].get("p_objective", 0.0) * mcj["doubled_sigma"]["D"]["n"])
    K["M2_N"] = "%d" % mcj["doubled_sigma"]["A"]["n"]
    for c in "ABCD":
        K["CAL_P" + c] = f(mcj["calibration"][c]["p_window"], 3)
        K["CAL_P05_" + c] = f(mcj["calibration"][c]["window_p05"], 2)
    K["CAL_N"] = "%d" % mcj["calibration"]["A"]["n"]
    cu = mcj["calibration_uncertainty"]
    K["CAL_SIG_KW"] = f(cu["sigma_k_wind"], 2)
    for c in "AD":
        K["B25_P" + c] = f(mcj["beta_25"][c]["p_window"], 3)
    K["B25_N"] = "%d" % mcj["beta_25"]["A"]["n"]
    K["P_MIN_EPI"] = f(min(mcj["model_form"]["A"]["p_window"], mcj["calibration"]["A"]["p_window"]), 2)

    # --- analyses added after the second internal review (2026-11-10) ---
    def ci2(x):
        lo, hi = x["p_window_ci95"]
        return "[%.2f, %.2f]" % (lo, hi)
    for c in "ABCD":
        K["M2%s_CI" % c] = ci2(mcj["doubled_sigma"][c])
        K["MF_CI" + c] = ci2(mcj["model_form"][c])
        K["CAL_CI" + c] = ci2(mcj["calibration"][c])
    K["M2A_P05"] = f(mcj["doubled_sigma"]["A"]["window_p05"], 2)
    K["MA_PLO"] = f(mcj["with_wind"]["A"]["p_window_ci95"][0], 3)
    lo, hi = mcj["with_wind"]["A"]["window_p05_ci95"]
    K["MA_P05_LO"], K["MA_P05_HI"] = f(lo, 2), f(hi, 2)
    for key, pre in (("msis", "MS_"), ("model_form_upper", "MFU_")):
        for c in "ABCD":
            x = mcj.get(key, {}).get(c)
            K[pre + "P" + c] = f(x["p_window"], 3) if x else "--"
            K[pre + "CI" + c] = ci2(x) if x else "--"
    mu = mcj.get("model_form_upper_sigma", {})
    K["MFU_PCT"] = f(mu.get("upper_pct", float("nan")), 1)
    for key, k2 in (("ign_att_1deg", "ATT1_PD"), ("ign_att_2deg", "ATT2_PD")):
        x = mcj.get(key, {}).get("D")
        K[k2] = f(x["p_window"], 3) if x else "--"
    if K.get("ATT1_PD") == K.get("ATT2_PD"):
        K["ATT_TXT"] = "a window probability of %s for D in both cases (%s without the error)" % (K["ATT1_PD"], f(mcj["with_wind"]["D"]["p_window"], 3))
    else:
        K["ATT_TXT"] = "window probabilities of %s and %s for D" % (K["ATT1_PD"], K["ATT2_PD"])
    ig = mcj.get("d_ignition_state", {})
    K["IGN_ALT"] = f(ig.get("alt_km", float("nan")), 1)
    K["IGN_M"] = f(ig.get("mach", float("nan")), 2)
    K["IGN_Q"] = "%.0f" % ig.get("q_kpa", float("nan"))
    K["IGN_ALPHA"] = f(ig.get("alpha_deg", float("nan")), 2)
    vv = load("vsb30_validation")
    fl = [x for x in vv["flights"] if not x["calibration"]]
    sig1 = math.sqrt(sum((100.0 * x["mc_sigma_apogee_km"] / x["apogee_actual"]) ** 2 for x in fl) / len(fl))
    K["V_SIG1"] = f(sig1, 1)
    K["V_SIGX"] = f(5.0 / sig1, 1)
    try:
        tp = load("thrust_pressure")
        rows = []
        def tprow(lab, isp, b1, cc, ho):
            return "%s & %.1f & %.2f & %.1f / %.3f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\" % (
                lab, isp, b1, ho["apogee_rms_pct"], ho["p_upper"], cc["A"]["window_s"], cc["B"]["window_s"], cc["C"]["window_s"],
                cc["D"]["window_s"], cc["A"]["burnout_mach"])
        for de in ("0.00", "0.40", "0.50"):
            x = tp["cases"][de]
            cc = x["concepts"]
            tag = de.replace("0.", "")
            K["TP_ISP" + tag] = f(x["isp_vac_s"], 0)
            K["TP_A" + tag], K["TP_D" + tag] = f(cc["A"]["window_s"], 1), f(cc["D"]["window_s"], 1)
            lab = "none (constant thrust)" if de == "0.00" else "%s~m" % de
            rows.append(tprow(lab, x["isp_vac_s"], x["s31_burn_s"], cc, x["held_out"]))
        b = tp["baseline"]
        K["TP_ISPB"] = f(tp["baseline_isp_s"], 0)
        rows.insert(2, tprow("%.2f~m (baseline)" % tp.get("baseline_exit_diameter_m", 0.45),
                             tp["baseline_isp_s"], tp["baseline_s31_burn_s"], b, tp["baseline_held_out"]))
        K["TP_ROWS"] = "\n".join(rows)
    except FileNotFoundError:
        pass
    try:
        rl = load("robust_local")
        rlg = rl["grid"]
        K["RL_P05_LO"] = f(min(x["tau_use_p05"] for x in rlg), 1)
        K["RL_P05_HI"] = f(max(x["tau_use_p05"] for x in rlg), 1)
        bst = rl["best"]
        K["RL_BEST"] = "%.0f\\textdegree\\ and %.1f~s" % (bst["el"], bst["coast"])
        K["RL_N6"] = "%d" % sum(x["tau_use_p05"] >= 6.0 for x in rlg)
        K["RL_NG"] = "%d" % len(rlg)
        K["RL_BEST_C"], K["RL_SEL_C"] = f(rl["best_confirm_p05"], 1), f(rl["selected_confirm_p05"], 1)
        K["RL_DP05"] = sgn(rl["best_minus_selected_p05"], 2)
        lo, hi = rl["best_minus_selected_p05_ci95"]
        K["RL_DP05_LO"], K["RL_DP05_HI"] = sgn(lo, 2), sgn(hi, 2)
        K["RL_EL_LO"] = "%.0f" % min(x["x"]["el"] for x in rlg)
        K["RL_EL_HI"] = "%.0f" % max(x["x"]["el"] for x in rlg)
        K["RL_CO_LO"] = f(min(x["x"]["coast"] for x in rlg), 1)
        K["RL_CO_HI"] = f(max(x["x"]["coast"] for x in rlg), 1)
        K["RL_EDGE"] = "on the edge of" if rl["best_on_edge"] else "inside"
        if bst == rl["selected"]:
            K["RL_CONCL"] = "(the best grid point is the selected setting)."
        elif lo <= 0.0 <= hi:
            K["RL_CONCL"] = ("an interval that contains zero: the 5th percentile stays within a few tenths of a second along a ridge "
                             "on which a lower elevation calls for a shorter coast, so the selection identifies a robust ridge rather than a sharp optimum.")
        elif lo > 0.0:
            K["RL_CONCL"] = "a significant improvement, so the screened candidate set was too coarse and the local optimum is preferable."
        else:
            K["RL_CONCL"] = "a significant loss: the screened best point does not hold on independent draws (winner's bias)."
        # elevation x coast matrix: each cell "5th percentile of tau_use / P" (FT-18); selected setting in bold,
        # best grid point underlined
        els_ = sorted({x["x"]["el"] for x in rlg})
        cos_ = sorted({x["x"]["coast"] for x in rlg})
        cell = {(x["x"]["el"], x["x"]["coast"]): x for x in rlg}
        sel_ = rl["selected"]
        K["RL_NCOL"] = "c" * len(cos_)
        K["RL_HEAD"] = " & ".join("%.0f" % c_ for c_ in cos_)
        rows_ = []
        for e_ in els_:
            cs__ = []
            for c_ in cos_:
                x = cell.get((e_, c_))
                if x is None:
                    cs__.append("--")
                    continue
                t = "%.1f / %.2f" % (x["tau_use_p05"], x["p_req"])
                if abs(e_ - sel_["el"]) < 1e-9 and abs(c_ - sel_["coast"]) < 1e-9:
                    t = "\\textbf{%s}" % t
                if abs(e_ - bst["el"]) < 1e-9 and abs(c_ - bst["coast"]) < 1e-9:
                    t = "\\underline{%s}" % t
                cs__.append(t)
            rows_.append("%.0f\\textdegree & %s \\\\" % (e_, " & ".join(cs__)))
        K["RL_ROWS"] = "\n".join(rows_)
    except FileNotFoundError:
        pass
    try:
        mm = load("mass_ceiling_mc")
        K["MCM_N"] = "%d" % mm["n"]
        K["MCM_NBOOT"] = "%d" % mm["n_boot"]
        vals = []
        rows = []

        def ctext(r, grid):
            """Ceiling rounded inward to 10 kg; 'only a--b' when the band is below P_MIN at 400 kg
            (a grid-edge bound carries the sign of the edge)."""
            g0, g1 = grid[0], grid[-1]
            if r["ceiling_kg"] is not None:
                v = 10 * math.floor(r["ceiling_kg"] / 10.0)
                return ("$\\geq$%.0f" % v) if r["ceiling_kg"] >= g1 - 0.5 else "%.0f" % v
            if r["band_kg"]:
                a, b = r["band_kg"]
                sa = ("$\\leq$%.0f" % a) if a <= g0 + 0.5 else "%.0f" % (10 * math.ceil(a / 10.0))
                sb = ("$\\geq$%.0f" % b) if b >= g1 - 0.5 else "%.0f" % (10 * math.floor(b / 10.0))
                return "only %s--%s" % (sa, sb)
            return "none"
        for c in "ABCD":
            x = mm["concepts"][c]
            grid = [r["mass_kg"] for r in x["grid"]]
            for met, tag in (("window", ""), ("use", "U")):
                K["MCM_" + tag + c] = ctext(x[met], grid)
                if tag == "" and x[met]["ceiling_kg"] is not None:
                    vals.append(10 * math.floor(x[met]["ceiling_kg"] / 10.0))
            qd = x["use"].get("quadratic", {})
            qw = x["window"].get("quadratic", {})
            gof = ("%.2f / %.2f" % (qw["p_gof"], qd["p_gof"])) if qw and qd else "--"
            rows.append("%s & %.0f--%.0f & %.2f / %.2f & %.2f / %.2f & %s / %s & %s \\\\" % (
                c, grid[0], grid[-1], x["window"]["p_at_400"], x["use"]["p_at_400"], x["window"]["band_at_400"],
                x["use"]["band_at_400"], K["MCM_" + c], K["MCM_U" + c], gof))
        K["MCM_ROWS"] = "\n".join(rows)

        def _grp(cs):
            return cs[0] if len(cs) == 1 else ", ".join(cs[:-1]) + " and " + cs[-1]

        def _mcm_txt(tag):
            full, only = {}, {}
            for c in "ABCD":
                v = K["MCM_" + tag + c]
                (only if v.startswith("only") else full).setdefault(v.replace("only ", ""), []).append(c)
            a = " and ".join("up to %s~kg for %s" % (v, _grp(cs)) for v, cs in full.items())
            b = " and ".join("%s~kg for %s" % (v, _grp(cs)) for v, cs in only.items())
            return a, b, only
        wa, wb, _ = _mcm_txt("")
        ua, ub, uonly = _mcm_txt("U")
        t = "At the nominal settings, the lower 95\\% band of the probability of meeting the window stays at 0.90 or more from 400~kg " + wa
        if wb:
            t += ", whereas it reaches 0.90 only at " + wb
        t += "; with the link, it "
        t += ("reaches 0.90 only at " + ub) if ub else ""
        if any(c in cs for cs in uonly.values() for c in "BC"):
            t += ", because a heavier vehicle flies a lower trajectory whose window is closer to CLA"
        if ua:
            t += (", and stays at 0.90 or more " if ub else "stays at 0.90 or more from 400~kg ") + ua
        K["MCM_TXT"] = t
        # stage-mass sensitivity (run_stage_mass_sensitivity.py)
        sms = load("stage_mass_sensitivity")["cases"]
        base = [c for c in sms if c["case"] == "base"][0]
        alt = [c for c in sms if c["case"] != "base"]
        rms = [c["heldout_apogee_rms_pct"] for c in alt]
        K["SMS_RMS_LO"], K["SMS_RMS_HI"] = f(min(rms)), f(max(rms))
        isp = [c["isp"] for c in alt]
        K["SMS_ISP_LO"], K["SMS_ISP_HI"] = f(min(isp)), f(max(isp))
        dmax = {k: max(abs(c["concepts"][k]["window_s"] - base["concepts"][k]["window_s"]) for c in alt) for k in "ABCD"}
        K["SMS_DA"], K["SMS_DBC"], K["SMS_DD"] = f(dmax["A"]), f(max(dmax["B"], dmax["C"])), f(dmax["D"])
        _edge = any("leq" in K["MCM_" + t + c] or "geq" in K["MCM_" + t + c] for t in ("", "U") for c in "ABCD")
        K["MCM_EDGE"] = "; $\\leq$, $\\geq$: edge of the 20~kg mass grid" if _edge else ""
        K["MCM_RANGE"] = ("%.0f--%.0f" % (min(vals), max(vals))) if vals else "--"
    except FileNotFoundError:
        pass
    try:
        qb = load("q_bias")
        K["HQ_DB"] = f(qb["max_abs_change_s"], 2)
    except FileNotFoundError:
        pass
    # third-round additions
    K["B_APOR"] = "%.0f" % N["B"]["apogee_km"]
    _e = sorted(op["max_dwell_400"][c]["x"]["el"] for c in "BC")
    K["OP4_ELBC"] = ("%.0f" % _e[0]) if round(_e[0]) == round(_e[1]) else "%.0f--%.0f" % (_e[0], _e[1])
    K["B_WDR"] = "%.0f" % (10 * round(sum(N["B"]["window_downrange_km"]) / 20.0))
    K["B_WM"] = f(sum(N["B"]["window_mach"]) / 2.0, 1)
    K["A_SEPQ"] = "%.0f" % N["A"]["burnout_q_kpa"]
    K["D_BOH"] = f(N["D"]["burnout_alt_km"], 0)
    K["A_SEPM"] = f(N["A"]["burnout_mach"], 1)
    K["D_UNP"] = f(N["D"]["window_s"] - (N["D"]["t_burnout"] - N["D"]["window_t"][0]), 1)
    try:
        le = load("le_transient")["concepts"]
        cc = [le[c] for c in "ABCD" if le[c]]
        K["LE_TAU_LO"], K["LE_TAU_HI"] = f(min(x["f1.0"]["tau_s"] for x in cc), 1), f(max(x["f1.0"]["tau_s"] for x in cc), 1)
        K["LE_R1_LO"], K["LE_R1_HI"] = "%.0f" % (100 * min(x["f1.0"]["ratio_rise"] for x in cc)), "%.0f" % (100 * max(x["f1.0"]["ratio_rise"] for x in cc))
        K["LE_R3_LO"], K["LE_R3_HI"] = "%.0f" % (100 * min(x["f0.3"]["ratio_rise"] for x in cc)), "%.0f" % (100 * max(x["f0.3"]["ratio_rise"] for x in cc))
        K["LE_TE_LO"], K["LE_TE_HI"] = "%.0f" % min(x["f1.0"]["T_max_window_K"] for x in cc), "%.0f" % max(x["f1.0"]["T_max_window_K"] for x in cc)
        K["LE_TE3_LO"], K["LE_TE3_HI"] = "%.0f" % min(x["f0.3"]["T_max_window_K"] for x in cc), "%.0f" % max(x["f0.3"]["T_max_window_K"] for x in cc)
        K["LE_TAU3_HI"] = f(max(x["f0.3"]["tau_s"] for x in cc), 1)
        sol = [le[c]["f1.0"]["t_solidus_from_window_start_s"] for c in "ABCD" if le[c]]
        sol_ok = [x for x in sol if x is not None]
        K["LE_NSOL"] = ["none", "one", "two", "three", "all"][len(sol_ok)]
        K["LE_TSOL_LO"] = f(min(sol_ok), 1) if sol_ok else "--"
        K["LE_TSOL_HI"] = f(max(sol_ok), 1) if sol_ok else "--"
        for c in "ABCD":
            K["LE_QC_" + c] = "%.0f" % le[c]["f1.0"]["Q_cold_MJ"]
            K["LE_QT_" + c] = "%.0f" % le[c]["f1.0"]["Q_transient_MJ"]
        K["LE_QT_MAX"] = "%.0f" % max(le[c]["f1.0"]["Q_transient_MJ"] for c in "ABCD")
        K["LE_ROWS"] = "\n".join("%s & %.1f & %.0f & %.0f & %.0f & %s & %.0f / %.0f & %.0f & %.0f \\\\" % (
            c, le[c]["f1.0"]["tau_s"], le[c]["f1.0"]["T_max_window_K"], le[c]["f1.0"]["Teq_max_window_K"],
            100 * le[c]["f1.0"]["ratio_rise"],
            f(le[c]["f1.0"]["t_solidus_from_window_start_s"], 1) if le[c]["f1.0"]["t_solidus_from_window_start_s"] is not None else "--",
            le[c]["f1.0"]["Q_cold_MJ"], le[c]["f1.0"]["Q_transient_MJ"], le[c]["f0.3"]["T_max_window_K"],
            100 * le[c]["f0.3"]["ratio_rise"]) for c in "ABCD" if le[c])
    except FileNotFoundError:
        pass
    try:
        ig = load("ignition_check")["cases"]
        z = [x for x in ig if x["angle_deg"] == 0.0 and x["cn_scale"] == 1.0][0]
        k5 = [x for x in ig if x["angle_deg"] == 5.0 and x["cn_scale"] == 1.0]
        k5h = [x for x in ig if x["angle_deg"] == 5.0 and x["cn_scale"] == 0.5]
        zh = [x for x in ig if x["angle_deg"] == 0.0 and x["cn_scale"] == 0.5][0]
        K["IG_T5"] = f(max(x["t_alpha_below_0p5_s"] for x in k5), 1)
        K["IG_T5H"] = f(max(x["t_alpha_below_0p5_s"] for x in k5h), 1)
        K["IG_DW"] = f(max(max(abs(x["window_s"] - z["window_s"]) for x in ig if x["cn_scale"] == 1.0),
                           max(abs(x["window_s"] - zh["window_s"]) for x in ig if x["cn_scale"] == 0.5)), 1)
        K["IG_D0H"] = f(zh["window_s"], 1)
        K["IG_QA5"] = "%.0f" % max(x["qalpha_max_kpa_deg"] for x in k5 + k5h)
        K["IG_ROWS"] = "\n".join("%.0f\\textdegree & %s & %.2f & %.1f & %s & %.0f & %.2f & %.2f \\\\" % (
            x["angle_deg"], "%.0f\\textdegree" % (x["phase"] * 180 / math.pi) if x["angle_deg"] > 0 else "--",
            x["cn_scale"], x["alpha_max_deg"], f(x["t_alpha_below_0p5_s"], 1) if x["angle_deg"] > 0 else "--",
            x["qalpha_max_kpa_deg"], x["burnout_mach"], x["window_s"]) for x in ig)
    except FileNotFoundError:
        pass
    try:
        to = load("tipoff_check")
        K["TO_SIG"] = f(to["sigma_tip_deg"], 2)
        K["TO_N"] = "%d" % to["N"]
        for c in "ABD":
            x = to["concepts"][c]
            K["TO_PU_" + c] = f(x["p_use"], 3)
            K["TO_PW_" + c] = f(x["p_window"], 3)
            K["TO_R3_" + c] = "%.0f" % x["r3_km"]
        K["TO_ROWS"] = "\n".join("%s & %.3f [%.3f, %.3f] & %.3f [%.3f, %.3f] & %.2f & %.0f \\\\" % (
            c, to["concepts"][c]["p_window"], *to["concepts"][c]["p_window_ci95"], to["concepts"][c]["p_use"],
            *to["concepts"][c]["p_use_ci95"], to["concepts"][c]["tau_use_p05"], to["concepts"][c]["r3_km"]) for c in "ABD")
    except FileNotFoundError:
        pass
    try:
        sb = load("site_bias")
        K["SB_CLA"], K["SB_ESR"] = sgn(sb["CLA"]["mean_err_msis_pct"], 1), sgn(sb["Esrange"]["mean_err_msis_pct"], 1)
        K["SB_RMS"] = f(sb["rms_msis_pct"], 1)
        K["SB_DV01"] = f(sb["delta_v01_km"], 1)
        K["SB_ROWS"] = "\n".join("%s & %s & %.1f & %.1f & %s & %s & %s \\\\" % (
            x["flight"], x["site"], x["actual_km"], x["us76_km"], sgn(x["delta_km"], 1), sgn(x["err_us76_pct"], 1),
            sgn(x["err_msis_pct"], 1)) for x in sb["flights"])
    except FileNotFoundError:
        pass
    try:
        tp = load("thrust_pressure")
        aa = [tp["cases"][d]["concepts"]["A"]["window_s"] for d in ("0.40", "0.50")]
        dd = [tp["cases"][d]["concepts"]["D"]["window_s"] for d in ("0.40", "0.50")]
        rng = lambda v: f(min(v), 1) if f(min(v), 1) == f(max(v), 1) else "%s--%s" % (f(min(v), 1), f(max(v), 1))  # noqa: E731
        K["TP_ARNG"], K["TP_DRNG"] = rng(aa), rng(dd)
    except FileNotFoundError:
        pass
    # ---- fourth-round additions ----
    from shyma import vehicle as _veh
    _A = math.pi * 0.557 ** 2 / 4.0
    K["BETA_TV"] = "%.0f" % (10 * round(400.0 / (_A * _veh.cd_generic(6.0)) / 10.0))
    K["D_BETA_STACK"] = "%.0f" % (100 * round((400.0 + _veh.S30_DRY) / (_A * _veh.cd_generic(6.0)) / 100.0))
    K["A_SEPM2"] = f(N["A"]["burnout_mach"], 2)
    ej = N["D"].get("expose_state", {})
    K["D_JET_H"], K["D_JET_M"] = f(ej.get("alt_km", float("nan")), 1), f(ej.get("mach", float("nan")), 1)
    K["D_JET_Q"], K["D_JET_G"] = "%.0f" % ej.get("q_kpa", float("nan")), f(ej.get("thrust_accel_g", float("nan")), 1)
    # NESC 17
    m17 = nesc["case17"]["metrics"]
    K["N17X"] = f(max(m17[q]["max_abs_err"] / m17[q]["max_spread"] for q in ("alt_m", "ve_mps", "vd_mps")), 1)
    ins = nesc["case17"].get("init_sensitivity", {})
    for nm, kk in (("SIM 04", "N17X04"), ("SIM 05", "N17X05"), ("SIM 06", "N17X06")):
        if nm in ins:
            K[kk] = f(ins[nm]["alt_m"]["max_abs_err"] / ins[nm]["alt_m"]["max_spread"], 1)
    # JSBSim powered
    try:
        jt = load("jsbsim_thrust")
        K["JT_MACH"] = f(100 * jt["max_rel_mach"], 3)
        K["JT_Q"] = f(100 * jt["max_rel_q"], 3)
        K["JT_MAX"] = f(100 * max(jt["max_rel_mach"], jt["max_rel_q"]), 3)
        K["JT_DT"] = f(1000 * jt["max_dt_s"], 1)
    except FileNotFoundError:
        pass
    # residual tests
    zz = load("vsb30_zstats")
    K["Z_P2"] = f(zz["p_two_sided"], 2)
    K["Z_PMEAN"] = f(zz["p_mean"], 2)
    K["Z_MEAN"] = sgn(zz["z_mean"], 2)
    pw = zz["power_upper_chi2_by_sigma_d_fraction"]
    K["Z_POW05"], K["Z_POW1"], K["Z_POW15"] = f(pw["0.5"], 2), f(pw["1.0"], 2), f(pw["1.5"], 2)
    K["SD_UP"] = f(zz["sigma_d_upper95_pct"], 1)
    K["SD_UPC"] = f(zz["sigma_d_upper95_with_cal_pct"], 1)
    K["ZS_CAL"] = f(zz["sigma_cal_apogee_km"])
    vvj = load("vsb30_validation")
    flv = [x for x in vvj["flights"] if not x["calibration"]]
    s1 = math.sqrt(sum((100.0 * x["mc_sigma_apogee_km"] / x["apogee_actual"]) ** 2 for x in flv) / len(flv))
    sc = 100.0 * zz["sigma_cal_apogee_km"] / (sum(x["apogee_actual"] for x in flv) / len(flv))
    K["V_RMS_PERFECT"] = f(math.sqrt(s1 ** 2 + sc ** 2), 1)
    sbj = load("site_bias")
    K["SB_PLOW"] = f(sbj["p_lower_msis"], 3)
    K["SB_CHI2"] = f(sbj["chi2_msis"], 1)
    # thrust model variants: held-out statistics
    tpj = load("thrust_pressure")
    h0 = tpj["cases"]["0.00"]["held_out"]
    K["TP_HO_P00"], K["TP_HO_RMS00"] = f(h0["p_upper"], 3), f(h0["apogee_rms_pct"], 1)
    # multi-flight calibration and leave-one-out spread
    import statistics as _st
    lo_ = vvj["loo"]
    K["LOO_SISP"], K["LOO_SB1"] = f(_st.stdev([x["isp"] for x in lo_]), 1), f(_st.stdev([x["b1"] for x in lo_]), 1)
    try:
        mf = load("multiflight_calibration")
        K["MFC_ISP"], K["MFC_SISP"] = f(mf["isp"], 1), f(mf["sigma_isp"], 1)
        K["MFC_B1"], K["MFC_SB1"] = f(mf["b1"], 2), f(mf["sigma_b1"], 2)
        K["MFC_RHO"] = f(mf["corr"], 2)
        K["MFC_RMS"] = f(mf["apogee_rms_pct"], 1)
        K["MFC_LOO"] = f(mf["loo_apogee_rms_pct"], 1)
        K["MFC_RCHI2"] = f(mf["chi2"] / mf["dof"], 2)
        K["MFC_LSISP"], K["MFC_LSB1"] = f(mf["loo_sigma_isp"], 1), f(mf["loo_sigma_b1"], 2)
        K["MFC_DW"] = ", ".join(f(mf["concepts"]["us76"][c]["window_s"], 1) for c in "ABCD")
        K["MFC_DWM"] = ", ".join(f(mf["concepts"]["msis_cla"][c]["window_s"], 1) for c in "ABCD")
    except FileNotFoundError:
        pass
    K["CAL_PRHO"] = f(mcj["calibration_uncertainty"]["corr"], 2)
    # HyShot from the apogee
    ap = hy["sensitivity"].get("apogee_start")
    if ap:
        K["HY_AP_M35"] = f(ap["mach_35"], 2)
        K["HY_AP_EM35"] = f(100 * (ap["mach_35"] - 7.75) / 7.75, 1)
        K["HY_AP_EQ35"] = "%.0f" % (100 * (ap["q_35_kpa"] - 21.0) / 21.0)
    # Monte Carlo: useful test time with the link, epistemic spread, paired tests
    ww = mcj["with_wind"]
    for c in "ABCD":
        x = ww[c]
        lo, hi = x["p_use_ci95"]
        K["MU_P" + c] = f(x["p_use"], 3)
        K["MU2_P" + c] = f(x["p_use"], 2)
        K["MU_CI" + c] = "[%s, %s]" % (f(lo, 3), f(hi, 3))
        K["MU_PCI" + c] = "%s %s" % (K["MU_P" + c], K["MU_CI" + c])
        pass
    _lo, _hi = f(min(ww["B"]["p_use"], ww["C"]["p_use"]), 2), f(max(ww["B"]["p_use"], ww["C"]["p_use"]), 2)
    K["MU_PBC"] = _lo if _lo == _hi else "%s--%s" % (_lo, _hi)
    K["MU_LOSS_BC"] = "%.0f" % (100 * max(ww[c]["p_window"] - ww[c]["p_use"] for c in "BC"))   # % of runs losing the link
    for c in "ABCD":
        x = ww[c]
        K["M%s_PPX" % c] = f(x["p_window"], 3)
        K["MU%s_P05" % c], K["MU%s_P50" % c] = f(x["tau_use_p05"], 2), f(x["tau_use_p50"], 2)
    EPI = [("doubled_sigma", r"doubled $\sigma$"), ("model_form", "model discrepancy"),
           ("model_form_upper", "discrepancy, upper bound"), ("calibration", "calibration, V01"),
           ("calibration_loo", "calibration, single-flight spread"),
           ("calibration_multiflight", "calibration, ten flights$^{\\S}$"),
           ("msis", "NRLMSISE-00 (CLA)")]
    pv = mcj.get("paired_vs_base", {})
    rows = []
    for key, lab in EPI:
        if key not in mcj:
            continue
        cells = []
        for c in "ABCD":
            x = mcj[key][c]
            star = r"$^{\ast}$" if pv.get(key, {}).get(c, {}).get("use_ok", {}).get("p_mcnemar", 1.0) < 0.05 else ""
            cells.append(f(x["p_use"], 3) + star)
        rows.append(r"\quad %s & %s \\" % (lab, " & ".join(cells)))
    K["MC_EPI_ROWS"] = "\n".join(rows)
    for c in "ABCD":
        # range over the epistemic alternatives only (the base set is reported separately)
        vals = [mcj[k][c]["p_use"] for k, _ in EPI if k in mcj]
        K["EPI_" + c] = "%s--%s" % (f(min(vals), 2), f(max(vals), 2))
        lowk = min((k for k, _ in EPI if k in mcj), key=lambda k: mcj[k][c]["p_use"])
        K["EPI_LOWSET_" + c] = dict(EPI)[lowk]
    K["EPI_NALT"] = "%d" % sum(1 for k, _ in EPI if k in mcj)
    if "calibration_loo_uncertainty" in mcj:   # share of clipped t_b1 draws (4-20 s) in the single-flight-spread set
        from scipy.stats import norm
        cu_ = mcj["calibration_loo_uncertainty"]
        K["CL_CLIP"] = f(100 * (norm.cdf((4.0 - cu_["b1"]) / cu_["sigma_b1"]) + norm.sf((20.0 - cu_["b1"]) / cu_["sigma_b1"])), 1)
    # selection: primary criterion and epistemic re-flights
    rsj = load("robust_selection")
    _bc = rsj["concepts"]["B"]["confirm"]
    K["RS_R_B84"] = "%.0f" % _bc["nominal"]["range_mean_km"]
    K["RS_R_B88"] = "%.0f" % _bc["R2"]["range_mean_km"] if "R2" in _bc and _bc["R2"] else "--"
    for c_ in "ABCD":
        cp_ = rsj["concepts"][c_]["confirm"]["primary"]
        cn_ = rsj["concepts"][c_]["confirm"]["nominal"]
        K["RSP_EL_" + c_] = ("%.0f\\textdegree{}/%.0f~s" % (cp_["x"]["el"], cp_["x"]["coast"])) if c_ == "D" else "%.0f\\textdegree{}" % cp_["x"]["el"]
        K["RSP_P05_" + c_] = f(cp_["tau_use_p05"], 1)
        K["RSP_P_" + c_] = f(cp_["p_req"], 3)
        from scipy.stats import beta as _beta
        _n = int(rsj["n_confirm"]); _k = int(round(cp_["p_req"] * _n))
        K["RSP_LB_" + c_] = f(_beta.ppf(0.025, _k, _n - _k + 1) if _k > 0 else 0.0, 3)
        K["RS_NCONF"] = "%d" % _n
        K["RSP2_P_" + c_] = f(cp_["p_req"], 2)
        K["RSN_P05_" + c_] = f(cn_["tau_use_p05"], 1)
        K["RSN_P_" + c_] = f(cn_["p_req"], 3)
        K["RSP_R_" + c_] = "%.0f" % cp_["range_mean_km"]
    ap_ = rsj["concepts"]["A"]
    K["RS_AP_EL"] = "%.0f" % ap_["primary_x"]["el"]
    K["RS_AP_P"] = f(ap_["confirm"]["primary"]["p_req"], 3)
    dd = rsj["concepts"]["D"]
    dp = dd["confirm"]["primary"]
    K["RS_DP_EL"], K["RS_DP_CO"] = "%.0f" % dp["x"]["el"], "%.0f" % dp["x"]["coast"]
    K["RS_DP_P05"], K["RS_DP_P50"], K["RS_DP_R"] = f(dp["tau_use_p05"], 1), f(dp["tau_use_p50"], 1), "%.0f" % dp["range_mean_km"]
    nf = round((1.0 - dp["p_req"]) * dp["n"])
    K["RS_DP_FAIL"] = "no failure" if nf == 0 else ("%d failure%s" % (nf, "" if nf == 1 else "s"))
    K["RS_DP_PLO"] = f(dp["p_req_ci95"][0], 3)
    K["RS_D1_EL"], K["RS_D1_CO"] = "%.0f" % dd["R1_x"]["el"], "%.0f" % dd["R1_x"]["coast"]
    if dd["R1_x"] == dd["primary_x"] and dd["R3_x"] == dd["primary_x"]:
        K["RS_DSINGLE"] = "(the same setting maximizes both the lower bound of $P$ and the 5th percentile)"
    else:
        K["RS_DSINGLE"] = "(the single criteria give %.0f\\textdegree/%.0f~s and %.0f\\textdegree/%.0f~s)" % (
            dd["R1_x"]["el"], dd["R1_x"]["coast"], dd["R3_x"]["el"], dd["R3_x"]["coast"])
    try:
        re_ = load("robust_epistemic")
        K["RE_N"] = "%d" % re_["n"]
        st_ = re_["settings"]
        keyof = lambda x: "%.0f/%.0f" % (x["el"], x["coast"])  # noqa: E731
        nomk, prk = keyof(st_["nominal"]), keyof(st_["primary"])
        b = re_["sets"]["base"]
        K["RE_DIFF_BASE"] = sgn(b["primary_minus_nominal_p05"], 1)
        K["RE_DIFF_BASE_LO"], K["RE_DIFF_BASE_HI"] = sgn(b["primary_minus_nominal_p05_ci95"][0], 1), sgn(b["primary_minus_nominal_p05_ci95"][1], 1)
        others = [k for k in re_["sets"] if k != "base"]
        K["RE_NSETS"] = "%d" % len(others)
        K["RE_NWIN"] = "%d" % sum(re_["sets"][k]["primary_minus_nominal_p05"] > 0 for k in others)
        K["RE_NPWIN"] = "%d" % sum(re_["sets"][k][prk]["p"] > re_["sets"][k][nomk]["p"] for k in others)
        _w = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight"}
        K["RE_NSETS_W"] = _w.get(len(others), "%d" % len(others))
        K["RE_PWIN_TXT"] = "in all of them" if K["RE_NPWIN"] == K["RE_NSETS"] else "in %s of them" % _w.get(int(K["RE_NPWIN"]), K["RE_NPWIN"])
        ps = [re_["sets"][k][prk]["p"] for k in others]
        pn = [re_["sets"][k][nomk]["p"] for k in others]
        K["RE_PSEL_BASE"], K["RE_PNOM_BASE"] = f(b[prk]["p"], 3), f(b[nomk]["p"], 3)
        lose = [k for k in others if re_["sets"][k]["primary_minus_nominal_p05"] <= 0]
        if not lose:
            K["RE_OTHERS"] = ""
        elif all(re_["sets"][k][prk]["tau_use_p05"] == 0 and re_["sets"][k][nomk]["tau_use_p05"] == 0 for k in lose):
            K["RE_OTHERS"] = "(in the others both are zero)"
        else:
            K["RE_OTHERS"] = "(in the others the selected setting is not worse)"
        p05s = [re_["sets"][k][prk]["tau_use_p05"] for k in others]
        K["RE_P05_EPI"] = "%s--%s" % (f(min(p05s), 1), f(max(p05s), 1))
        K["RE_PSEL_LO"], K["RE_PSEL_HI"] = f(min(ps), 3), f(max(ps), 3)
        K["RE_PNOM_LO"], K["RE_PNOM_HI"] = f(min(pn), 3), f(max(pn), 3)
        setk = [k for k in re_["settings"] if re_["settings"][k]]
        uniq = []
        for k in ("nominal", "primary", "R1", "R3"):
            x = st_.get(k)
            if x and keyof(x) not in [u for u, _ in uniq]:
                uniq.append((keyof(x), k))
        labs = {"base": "base", "calibration": "calibration, V01", "calibration_loo": "calibration, single-flight spread",
                "calibration_multiflight": "calibration, ten flights", "doubled_sigma": r"doubled $\sigma$",
                "msis": "NRLMSISE-00 (CLA)"}
        head = "Set & " + " & ".join("%s (%s)" % (u, {"nominal": "nominal", "primary": "primary", "R1": "max.\\ $P$", "R3": "max.\\ 5th"}[k]) for u, k in uniq) + r" & $\Delta$5th, primary$-$nominal \\"
        body = []
        for kset in re_["sets"]:
            r = re_["sets"][kset]
            cells = ["%s / %s" % (f(r[u]["p"], 3), f(r[u]["tau_use_p05"], 1)) for u, _ in uniq]
            body.append("%s & %s & %s [%s, %s] \\\\" % (labs.get(kset, kset), " & ".join(cells), sgn(r["primary_minus_nominal_p05"], 1),
                                                       sgn(r["primary_minus_nominal_p05_ci95"][0], 1), sgn(r["primary_minus_nominal_p05_ci95"][1], 1)))
        K["RE_TABLE"] = ("\\begin{tabular}{l" + "c" * (len(uniq) + 1) + "}\n\\toprule\n" + head + "\n\\midrule\n" +
                         "\n".join(body) + "\n\\bottomrule\n\\end{tabular}")
    except FileNotFoundError:
        pass
    K["GD_450"] = f(cs["D_heavy"]["450"]["window_s"], 1)
    K["GD_500"] = f(cs["D_heavy"]["500"]["window_s"], 2)
    K["GD_500TXT"] = ("no setting meets the window at 500~kg" if cs["D_heavy"]["500"]["window_s"] == 0.0
                      else "%s~s at 500~kg" % K["GD_500"])
    obj = [(P, p, el) for P in "ABC" for i, p in enumerate(g["payloads"])
           for j, el in enumerate(g["elevations"]) if g[P][i][j]["objective_met"]]
    K["GO_N"] = "%d" % len(obj)
    if obj:
        K["GO_PMAX"] = "%d" % max(p for _, p, _ in obj)
        K["GO_EL"] = "%d--%d\\textdegree" % (min(el for _, _, el in obj), max(el for _, _, el in obj))
        rr = [g[P][g["payloads"].index(p)][g["elevations"].index(el)]["impact_range_km"] for P, p, el in obj]
        K["GO_R"] = "%.0f--%.0f" % (min(rr), max(rr))
    if obj:
        ps = sorted(set(p for _, p, _ in obj))
        cc = sorted(set(P for P, _, _ in obj))
        K["GO_TXT"] = "Concept%s %s, test-vehicle masses of %d kg or less" % (
            "s" if len(cc) > 1 else "", " and ".join(cc), ps[-1])
    else:
        K["GO_TXT"] = "none of the concepts at any mass or elevation"

    # thrust-curve shape (scripts/run_thrust_shape.py)
    ts = load("thrust_shape")
    tv = {k: v for k, v in ts["variants"].items() if k != "constant"}
    K["TS_K"] = "%.0f" % (100 * ts["k"])
    for c in "AD":
        nom = [v[c + "_nominal"]["dwell_s"] for v in tv.values()]
        ret = [v[c + "_retuned"]["dwell_s"] for v in tv.values()]
        K["TS_%s_LO" % c], K["TS_%s_HI" % c] = f(min(nom), 1), f(max(nom), 1)
        K["TS_%sR_LO" % c], K["TS_%sR_HI" % c] = f(min(ret), 1), f(max(ret), 1)
    K["TS_ISP_LO"] = f(min(v["isp"] for v in tv.values()), 1)
    K["TS_ISP_HI"] = f(max(v["isp"] for v in tv.values()), 1)
    K["TS_BM_LO"] = f(min(v["A_nominal"]["burnout_mach"] for v in tv.values()), 2)
    K["TS_BM_HI"] = f(max(v["A_nominal"]["burnout_mach"] for v in tv.values()), 2)
    rows = []
    for name, v in ts["variants"].items():
        rows.append("%s & %s & %s & %s & %s & %s (%s\\textdegree, %s~s) \\\\" % (
            name, f(v["isp"], 1), f(v["b1"], 2), f(v["A_nominal"]["dwell_s"], 2),
            f(v["D_nominal"]["dwell_s"], 2), f(v["D_retuned"]["dwell_s"], 2),
            "%.0f" % v["D_retuned"]["el"], "%.0f" % v["D_retuned"]["coast"]))
    K["TS_ROWS"] = "\n".join(rows)
    # pseudo-5-DOF attitude verification (scripts/run_weathercock_check.py)
    wc = load("weathercock_check")
    K["WC_ERR"] = f(100 * wc["small"]["max_error_rel"], 2)
    K["WC_ERR5"] = f(100 * wc["large"]["max_error_rel"], 1)
    K["WC_PER"] = f(wc["small"]["period_s"], 2)
    # rank-correlation sensitivity of the Monte Carlo dwell (scripts/run_mc_sensitivity.py)
    ms = load("mc_sensitivity")
    lab = ms["labels"]
    rows = []
    for c in "ABCD":
        top = ms["concepts"][c]["ranking"][:4]
        rows.append("%s & %s \\\\" % (c, "; ".join("%s (%+.2f)" % (lab[k], ms["concepts"][c]["rho"][k]) for k in top)))
    K["MS_ROWS"] = "\n".join(rows)
    K["MS_RHOMAX"] = f(max(abs(r) for c in "ABCD" for r in ms["concepts"][c]["rho"].values()), 2)
    atm = ("dT0", "dT1", "rho")
    K["MS_ATM_A"] = "%.0f" % (100 * sum(ms["concepts"]["A"]["share"].get(k, 0) for k in atm))
    K["MS_ATM_D"] = "%.0f" % (100 * sum(ms["concepts"]["D"]["share"].get(k, 0) for k in atm))
    K["MS_SUM_A"] = f(sum(r * r for r in ms["concepts"]["A"]["rho"].values()), 2)
    K["MS_SUM_D"] = f(sum(r * r for r in ms["concepts"]["D"]["rho"].values()), 2)
    K["MS_DRAG_D"] = "%.0f" % (100 * ms["concepts"]["D"]["share"]["cd"])
    # Monte Carlo
    mc = load("montecarlo")
    K["MC_WALL"] = "%.0f" % ((mc["cpu_s"] + mc.get("epistemic_cpu_s", 0.0)) / 60.0)
    for P in "ABCD":
        w, nw, x2 = mc["with_wind"][P], mc["no_wind"][P], mc["doubled_sigma"][P]
        Q = "M" + P
        K[Q + "_PP"] = f(w["p_window"], 3)
        lo, hi = w["p_window_ci95"]
        K[Q + "_PCIB"] = "[%s, %s]" % (f(lo, 3), f(hi, 3))
        K[Q + "_PCI"] = "%s %s" % (f(w["p_window"], 3), K[Q + "_PCIB"])
        K[Q + "_PL"] = f(w["p_window_link"], 3)
        K[Q + "_P05"] = f(w["window_p05"], 2)
        K[Q + "_P05CI"] = "[%s, %s]" % (f(w["window_p05_ci95"][0], 1), f(w["window_p05_ci95"][1], 1))
        K[Q + "_P50"] = f(w["window_p50"], 2)
        K[Q + "_P95"] = f(w["window_p95"], 2)
        K[Q + "_APO"] = r"$%.1f \pm %.1f$" % (w["apogee_mean"], w["apogee_std"])
        K[Q + "_R"] = "%.0f" % w["range_mean"]
        e = w["ellipse"]
        K[Q + "_AX"] = "%.1f/%.1f" % (e["sigma_major_km"], e["sigma_minor_km"])
        K[Q + "_R3"] = "%.0f" % e["r3sig_gauss_km"]
        K[Q + "_R99"] = f(e["r99_km"])
        K[Q + "_QA99"] = "%.0f" % w["qalpha_p99"]
        K[Q + "_3SNW"] = "%.0f" % (3 * nw["ellipse"]["sigma_major_km"])
        K[Q + "_FAC"] = f(e["sigma_major_km"] / nw["ellipse"]["sigma_major_km"])
        K["M2%s_P" % P] = f(x2["p_window"], 3)

    K = {k: (v + "{}" if isinstance(v, str) and v.endswith("\\textdegree") else v) for k, v in K.items()}
    left = []
    for pdir in ("paper_acta",):          # Acta Astronautica version (the earlier IEEE Access version in paper/ is frozen)
      for stem in ("main", "supplement"):
        tp = os.path.join(ROOT, pdir, stem + ".tpl.tex")
        if not os.path.exists(tp):
            continue
        src = open(tp).read()
        out = re.sub(r"@@([A-Z0-9_]+)@@", lambda m: K.get(m.group(1), m.group(0)), src)
        left += sorted(set(re.findall(r"@@([A-Z0-9_]+)@@", out)))
        with open(os.path.join(ROOT, pdir, stem + ".tex"), "w") as fo:
            fo.write(out)
    with open(os.path.join(ROOT, "results", "manuscript_values.json"), "w") as fo:
        json.dump(K, fo, indent=1)
    for rp in (os.path.join(ROOT, "docs", "review_round1", "RESPONSE.tpl.md"),
               os.path.join(ROOT, "docs", "review_round2", "RESPONSE_ROUND2.tpl.md")):
      if os.path.exists(rp):
        txt = open(rp).read()
        plain = {k: re.sub(r"\\textdegree(\{\})?|\\", "", v).replace("$", "").replace("--", "–") for k, v in K.items()}
        txt = re.sub(r"@@([A-Z0-9_]+)@@", lambda m: plain.get(m.group(1), m.group(0)), txt)
        left += sorted(set(re.findall(r"@@([A-Z0-9_]+)@@", txt)))
        open(rp.replace(".tpl", ""), "w").write(txt)
    print("filled %d keys; unfilled: %s" % (len(K), left))
    return 1 if left else 0


if __name__ == "__main__":
    sys.exit(main())
