"""Case study: pre-Phase A trade space of a VSB-30-boosted scramjet flight
experiment launched from the Alcantara Launch Center (CLA).

Four mission concepts are screened against the same a-priori test window:
  A - ascent test: the test vehicle separates at S30 burnout and is tested
      during the unpowered ascent, after a delay DSEP for thrust tail-off,
      separation, shroud/inlet-cover jettison and settling;
  B - re-entry test (separated): the test vehicle separates at t = 60 s,
      is reoriented outside the atmosphere and tested nose-first on descent;
  C - re-entry test (captive): as B, but the spent S30 stays attached
      (HyShot / HIFiRE-1 practice); the stack is reoriented and tested;
  D - boosted test (depressed, delayed ignition): the test vehicle stays
      attached to the S30, whose ignition is delayed by a longer interstage
      coast on a depressed launch, so that the S30 accelerates the experiment
      through the window (HIFiRE-2 practice); booster thrust is allowed in the
      window (F_b >= 0), from S30 ignition until apogee.

Requirement windows:
  threshold  Mach 6-8, 23-35 km, q 20-90 kPa, >= 3 s, unpowered
  objective  Mach >= 6.75 (14-X B design point Mach 7), q 40-90 kPa, >= 3 s
  q-band dwell: time inside the threshold window with q in 40-60 kPa
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import analysis, earth, flight, vehicle  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
SITE = earth.STATIONS["CLA"]
AZIMUTH = 45.0
WINDOW = dict(mach=(6.0, 8.0), alt_km=(23.0, 35.0), q_kpa=(20.0, 90.0), min_duration=3.0)
OBJECTIVE = dict(mach=(6.75, 8.0), alt_km=(23.0, 35.0), q_kpa=(40.0, 90.0), min_duration=3.0)
STEADY = dict(mach=(6.0, 8.0), alt_km=(23.0, 35.0), q_kpa=(40.0, 60.0), min_duration=0.0)
RN_LE = 0.002     # leading-edge radius (m), unswept cylinder
RN_NOSE = 0.050   # nose-cap / shroud-nose radius (m), sphere
# Test-vehicle ballistic coefficient at Mach 6 (kg/m^2): a 400 kg vehicle of the VSB-30 diameter (0.557 m)
# with the generic slender drag curve (C_D = 0.22 at Mach 6), i.e. the drag of the generic payload nose,
# so that the attached stack and the separated vehicle are consistent (about 7460 kg/m^2).
BETA_TV = 400.0 / (math.pi * 0.557 ** 2 / 4.0 * vehicle.cd_generic(6.0))
DSEP = 2.0        # Concept A: burnout-to-test delay (s)
BETA_STAGE = 2500.0  # spent finned stages, nose-forward (C_D ~ 0.45 on 0.244 m^2)
BETA_TUMBLE = 300.0  # spent stages, tumbling (sensitivity)
NOMINAL = {"A": 63.0, "B": 84.0, "C": 84.0, "D": 72.5}
PHASE = {"A": "ascent", "B": "descent", "C": "descent", "D": "ascent"}
COAST_D = 24.0     # Concept D nominal interstage coast (s): optimum of the inverse design at 400 kg (run_optimizer.py)
COAST_STD = 2.0    # standard VSB-30 interstage coast (s)
CAPTIVE_SEP = 1.0e9


def t_burnout(v):
    return v.t_burnout


def build(concept, payload, beta_tv=BETA_TV, coast=None, **kw):
    if concept == "C":
        return vehicle.vsb30(payload, t_payload_sep=CAPTIVE_SEP, beta_tv=beta_tv, **kw)
    if concept == "D":
        return vehicle.vsb30(payload, t_payload_sep=CAPTIVE_SEP, beta_tv=beta_tv,
                             coast=COAST_D if coast is None else coast, **kw)
    v0 = vehicle.vsb30(payload, **{k: x for k, x in kw.items() if k in ("burn_factor", "s31_burn", "isp")})
    t_sep = v0.t_burnout if concept == "A" else 60.0
    return vehicle.vsb30(payload, t_payload_sep=t_sep, beta_tv=beta_tv, **kw)


def fly(concept, payload, el, az=AZIMUTH, rho_factor=1.0, dt=0.01, rec=2, wind=None,
        atm=None, dT=None, beta_tv=BETA_TV, ign_att=None, tip_att=None, **kw):
    """ign_att = (angle_deg, phase_rad): attitude error of the stack at S30 ignition;
    tip_att = (angle_deg, phase_rad): attitude error just after rail exit (tip-off)."""
    v = build(concept, payload, beta_tv=beta_tv, **kw)
    # end of 5-DOF ascent: burnout for A and D, t = 60 s (exo-atmospheric) for B and C
    t_sw = v.t_burnout if concept in "AD" else max(60.0, v.t_burnout + 1.0)
    kick = None if ign_att is None else (v.stages[1].t_ign, ign_att[0], ign_att[1])
    if tip_att is not None:
        kick = (None, tip_att[0], tip_att[1])
    tr = flight.fly(v, SITE, el, az, t_sw, wind=wind, atm=atm, dT=dT,
                    rho_factor=rho_factor, dt=dt, record_every=rec, att_kick=kick)
    return v, tr


def windows(concept, v, tr, dsep=DSEP):
    t_min = v.t_burnout + dsep if concept == "A" else None
    pw = concept == "D"
    if pw:
        t_min = [s for s in v.stages if s.name == "S30"][0].t_ign
    ph = PHASE[concept]
    w = analysis.TestWindow(phase=ph, t_min=t_min, powered=pw, **WINDOW).evaluate(tr)
    wo = analysis.TestWindow(phase=ph, t_min=t_min, powered=pw, **OBJECTIVE).evaluate(tr)
    ws = analysis.TestWindow(phase=ph, t_min=t_min, powered=pw, **STEADY).evaluate(tr)
    return w, wo, ws, t_min


def useful_window(tr, w, t_ready=0.0, need_link=True, station="CLA", ant_loss_db=10.0, fade_db=6.0,
                  mask_deg=3.0):
    """Useful test time: the longest continuous interval in which the
    flight conditions of the window w are met, the acquisition is ready
    (t_ready s after the window opens) and, if need_link, the station sees the
    vehicle above mask_deg with a degraded link margin >= 0 dB. With
    need_link=False (on-board recording with recovery) only the flight window
    and readiness count. Disjoint pieces are not added."""
    if w["t_start"] is None:
        return 0.0
    t0, t1 = w["t_start"] + t_ready, w["t_end"]
    if t1 <= t0:
        return 0.0
    if not need_link:
        return t1 - t0
    ks = [k for k in range(len(tr.t)) if t0 <= tr.t[k] <= t1]
    if not ks:
        return 0.0
    cov = analysis.coverage(tr, [station], idx=ks)[station]
    ok = [e >= mask_deg and m - ant_loss_db - fade_db >= 0.0 for e, m in zip(cov["el"], cov["margin"])]
    best, i = 0.0, 0
    while i < len(ks):
        if not ok[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(ks) and ok[j + 1]:
            j += 1
        a = t0 if i == 0 else tr.t[ks[i]]
        b = t1 if j == len(ks) - 1 else tr.t[ks[j]]
        best = max(best, b - a)
        i = j + 1
    return best


def summarize(concept, v, tr, full=False, dsep=DSEP, atm=None, light=False):
    w, wo, ws, t_min = windows(concept, v, tr, dsep)
    ta, ha, ia = tr.apogee()
    dr = analysis.downrange_km(tr, SITE)
    t_exp = t_min if concept == "A" else 60.0
    if concept == "D":   # shroud jettisoned 2 s before the test (HIFiRE-2 practice)
        t_exp = (w["t_start"] - 2.0) if w["t_start"] is not None else v.t_burnout
    if light:   # Monte Carlo: no heating post-processing
        qle = Qle = Tle = qhle = qn = [0.0] * len(tr.t)
    else:
        qle, Qle, Tle, qhle = analysis.heating_history(tr, RN_LE, geometry="cylinder", t_expose=t_exp, atm=atm)
        qn, Qn, Tn, qhn = analysis.heating_history(tr, RN_NOSE, geometry="sphere", atm=atm)
    out = {"window_s": w["duration"], "window_met": w["met"], "t_window": w["t_start"],
           "objective_s": wo["duration"], "objective_met": wo["met"],
           "steady_s": ws["duration"],
           "apogee_km": ha / 1e3, "t_apogee": ta, "impact_range_km": dr[-1],
           "impact_lat": tr.lat[-1], "impact_lon": tr.lon[-1], "t_impact": tr.t[-1],
           "q_max_kpa": max(tr.q) / 1e3, "mach_max": max(tr.mach)}
    asc = [k for k in range(len(tr.t)) if tr.t[k] <= ta]
    out["q_max_ascent_kpa"] = max(tr.q[k] for k in asc) / 1e3
    alpha = getattr(tr, "alpha", None)
    if alpha:
        out["qalpha_max_kpa_deg"] = max(tr.q[k] * alpha[k] for k in range(min(len(alpha), len(tr.t)))
                                        if tr.thrust[k] > 0 or tr.t[k] < 60.0) / 1e3
        out["alpha_max_ascent_deg"] = max(alpha[k] for k in range(min(len(alpha), len(tr.t)))
                                          if tr.t[k] > 1.0 and tr.q[k] > 5e3)
    # shroud / nose heating before exposure (boost phase), and burnout state
    pre = [k for k in range(len(tr.t)) if tr.t[k] <= (t_exp or 0.0)]
    out["qdot_shroud_nose_boost_MW"] = (max(qn[k] for k in pre) / 1e6) if pre else 0.0
    if concept in "AD" and t_exp is not None:   # state at separation / shroud jettison
        sj = analysis.window_state(tr, t_exp)
        acc = v.thrust(t_exp) / v.mass(t_exp) / 9.80665
        out["expose_state"] = {"t": t_exp, "alt_km": sj["alt"] / 1e3, "mach": sj["mach"], "q_kpa": sj["q"] / 1e3,
                               "thrust_accel_g": acc}
    st = analysis.window_state(tr, v.t_burnout)
    out.update({"t_burnout": v.t_burnout, "burnout_alt_km": st["alt"] / 1e3,
                "burnout_mach": st["mach"], "burnout_q_kpa": st["q"] / 1e3, "burnout_fpa": st["fpa"]})
    if w["t_start"] is not None:
        t0, t1 = w["t_start"], w["t_end"]
        s0, s1 = analysis.window_state(tr, t0), analysis.window_state(tr, t1)
        ks = [k for k in range(len(tr.t)) if t0 <= tr.t[k] <= t1]
        cov = analysis.coverage(tr, ["CLA", "CLBI"], idx=ks)
        out["window_stations"] = {}
        for n in ("CLA", "CLBI"):
            els = cov[n]["el"] or [float("nan")]
            mg = cov[n]["margin"] or [float("nan")]
            out["window_stations"][n] = {"el_min": min(els), "el_max": max(els),
                                         "margin_min_db": min(mg)}
        out["window_alt_km"] = [s0["alt"] / 1e3, s1["alt"] / 1e3]
        out["window_mach"] = [s0["mach"], s1["mach"]]
        out["window_q_kpa"] = [s0["q"] / 1e3, s1["q"] / 1e3]
        out["window_q_mean_kpa"] = sum(tr.q[k] for k in ks) / max(len(ks), 1) / 1e3
        out["window_fpa_entry_deg"] = s0["fpa"]
        out["window_downrange_km"] = [earth.great_circle_km(SITE[0], SITE[1], s0["lat"], s0["lon"]),
                                      earth.great_circle_km(SITE[0], SITE[1], s1["lat"], s1["lon"])]
        kw = ks or [0]
        out["qdot_le_window_MW"] = max(qle[k] for k in kw) / 1e6
        out["qdot_le_hw_window_MW"] = max(qhle[k] for k in kw) / 1e6
        out["T_le_window_K"] = max(Tle[k] for k in kw)
        out["Q_le_to_window_end_MJ"] = Qle[kw[-1]] / 1e6
        out["T0_window_K"] = max(
            analysis.T_from_h(analysis.h_air((tr.v[k] / tr.mach[k]) ** 2 / (1.4 * 287.05287))
                              + 0.5 * tr.v[k] ** 2)
            for k in kw)
    if full:
        cov = analysis.coverage(tr, ["CLA", "CLBI"])
        out["visible_fraction_flight"] = {n: analysis.visible_fraction(tr, cov, n) for n in ("CLA", "CLBI")}
        step = max(1, len(tr.t) // 1500)
        out["profile"] = {
            "t": tr.t[::step], "alt_km": [a / 1e3 for a in tr.alt[::step]],
            "mach": tr.mach[::step], "q_kpa": [q / 1e3 for q in tr.q[::step]],
            "downrange_km": dr[::step], "lat": tr.lat[::step], "lon": tr.lon[::step],
            "qdot_le_MW": [q / 1e6 for q in qle[::step]], "T_le_K": Tle[::step],
            "el_CLA": cov["CLA"]["el"][::step], "el_CLBI": cov["CLBI"]["el"][::step],
            "margin_CLA": cov["CLA"]["margin"][::step], "margin_CLBI": cov["CLBI"]["margin"][::step],
        }
        out["window_t"] = [w["t_start"], w["t_end"]]
    return out


def bearing(a, b):
    la1, lo1, la2, lo2 = map(math.radians, (a[0], a[1], b[0], b[1]))
    y = math.sin(lo2 - lo1) * math.cos(la2)
    x = math.cos(la1) * math.sin(la2) - math.sin(la1) * math.cos(la2) * math.cos(lo2 - lo1)
    return math.degrees(math.atan2(y, x)) % 360.0


def stage_impacts(concept, payload, el, az=AZIMUTH):
    v = build(concept, payload)
    s31 = [s for s in v.stages if s.name == "S31"][0]
    s30 = [s for s in v.stages if s.name == "S30"][0]
    out = {}
    for name, st in (("S31", s31), ("S30", s30)):
        if concept in "CD" and name == "S30":
            continue
        la, lo, t, _ = flight.spent_stage_impact(v, SITE, el, az, st.t_drop, BETA_STAGE)
        la2, lo2, _, _ = flight.spent_stage_impact(v, SITE, el, az, st.t_drop, BETA_TUMBLE)
        out[name] = {"range_km": earth.great_circle_km(SITE[0], SITE[1], la, lo), "lat": la, "lon": lo,
                     "range_tumbling_km": earth.great_circle_km(SITE[0], SITE[1], la2, lo2)}
    return out


if __name__ == "__main__":
    t0 = time.time()
    res = {"window": WINDOW, "objective": OBJECTIVE, "steady": STEADY, "azimuth": AZIMUTH,
           "isp": vehicle.ISP_VSB30, "rn_le": RN_LE, "rn_nose": RN_NOSE, "beta_tv": BETA_TV,
           "dsep": DSEP, "beta_stage": BETA_STAGE, "earth_model": earth.MODEL}
    res["coast_D"] = COAST_D
    res["nominal"] = {}
    for c in "ABCD":
        v, tr = fly(c, 400.0, NOMINAL[c], rec=1)
        res["nominal"][c] = dict(summarize(c, v, tr, full=True), elevation=NOMINAL[c], payload=400.0)
        res["nominal"][c]["stages"] = stage_impacts(c, 400.0, NOMINAL[c])
        print(c, {k: x for k, x in res["nominal"][c].items() if k != "profile"}, flush=True)
    # standard elevation for Concept A (burnout above the window)
    v, tr = fly("A", 400.0, 84.0, rec=1)
    res["A_at_84"] = {k: x for k, x in summarize("A", v, tr).items() if k.startswith("burnout") or k == "window_s"}
    t_nom = time.time() - t0
    # sensitivities of the nominal concepts
    sens = {}
    for cst in (2.0, 10.0, 16.0, 20.0, 24.0, 27.0, 30.0):
        v, tr = fly("D", 400.0, NOMINAL["D"], coast=cst)
        s = summarize("D", v, tr)
        sens["coast_D_%d" % cst] = {k: s[k] for k in ("window_s", "objective_s", "steady_s", "impact_range_km",
                                                     "q_max_ascent_kpa", "burnout_alt_km", "burnout_mach")}
    for b in (2500.0, 5000.0, 10000.0, 15000.0):
        for c in "ABCD":
            v, tr = fly(c, 400.0, NOMINAL[c], beta_tv=b)
            s = summarize(c, v, tr)
            sens["beta_%s_%d" % (c, b)] = {k: s[k] for k in ("window_s", "objective_s", "steady_s", "impact_range_km")}
    v, tr = fly("A", 400.0, NOMINAL["A"])
    for d in (0.0, 1.0, 2.0, 3.0, 4.0):
        s = summarize("A", v, tr, dsep=d)
        sens["dsep_%.0f" % d] = {k: s[k] for k in ("window_s", "objective_s", "steady_s")}
    az_clbi = bearing(SITE, earth.STATIONS["CLBI"])
    for c in "BC":
        v, tr = fly(c, 400.0, NOMINAL[c], az=az_clbi, rec=1)
        s = summarize(c, v, tr)
        sens["toward_CLBI_" + c] = {"azimuth": az_clbi, "window_s": s["window_s"],
                                    "impact_range_km": s["impact_range_km"],
                                    "stations": s.get("window_stations")}
    from shyma.atmosphere import msis_atmosphere
    msis = msis_atmosphere(SITE[0], SITE[1], "2021-12-14T15:00", f107=100.0, f107a=100.0, ap=5.0)
    for c in "ABCD":
        v, tr = fly(c, 400.0, NOMINAL[c], rec=1, atm=msis)
        s = summarize(c, v, tr, atm=msis)
        sens["msis_" + c] = {k: s[k] for k in ("window_s", "objective_s", "apogee_km", "impact_range_km",
                                               "burnout_mach", "burnout_alt_km")}
    res["sensitivity"] = sens
    # trade grid
    els = [54 + 2 * i for i in range(18)]          # 54..88 deg
    pls = [250 + 50 * i for i in range(8)]          # 250..600 kg
    grid = {"elevations": els, "payloads": pls}
    for c in "ABC":
        grid[c] = []
        for p in pls:
            row = []
            for el in els:
                v, tr = fly(c, p, el, dt=0.05, rec=2)
                s = summarize(c, v, tr)
                row.append({k: s[k] for k in ("window_s", "window_met", "objective_s", "objective_met",
                                              "steady_s", "impact_range_km", "apogee_km", "burnout_mach",
                                              "mach_max", "q_max_ascent_kpa")})
            grid[c].append(row)
        print("grid", c, "done", flush=True)
    # Concept D grid: launch elevation x interstage coast, for three masses
    els_d = [56 + 2 * i for i in range(13)]        # 56..80 deg
    coasts = [2 * i for i in range(16)]            # 0..30 s
    gd = {"elevations": els_d, "coasts": coasts, "payloads": [300, 350, 400]}
    for p in gd["payloads"]:
        rows = []
        for cst in coasts:
            row = []
            for el in els_d:
                v, tr = fly("D", float(p), el, dt=0.05, rec=2, coast=max(cst, 0.5))
                s = summarize("D", v, tr, light=True)
                row.append({k: s[k] for k in ("window_s", "window_met", "objective_s", "objective_met",
                                              "steady_s", "impact_range_km", "q_max_ascent_kpa")})
            rows.append(row)
        gd[str(p)] = rows
        print("grid D", p, "done", flush=True)
    grid["D"] = gd
    # heavier test vehicles for Concept D: best dwell over elevation x coast
    heavy = {}
    for p in (450, 500):
        best = (0.0, None, None)
        for el in range(60, 80, 2):
            for cst in range(0, 30, 2):
                v, tr = fly("D", float(p), el, dt=0.05, rec=2, coast=max(cst, 0.5))
                w = windows("D", v, tr)[0]["duration"]
                if w > best[0]:
                    best = (w, el, cst)
        heavy[str(p)] = {"window_s": best[0], "elevation": best[1], "coast": best[2]}
    res["D_heavy"] = heavy
    res["grid"] = grid
    res["cpu_s"] = {"nominal_flights": t_nom, "total": time.time() - t0,
                    "n_trajectories": 3 * len(els) * len(pls),
                    "n_trajectories_D": len(els_d) * len(coasts) * 3}
    with open(os.path.join(OUT, "casestudy.json"), "w") as f:
        json.dump(res, f)
    print(res["cpu_s"])
