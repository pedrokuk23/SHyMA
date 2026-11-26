"""Browser driver of the SHyMA reviewer workbench (runs under Pyodide).

It calls the unmodified SHyMA package and case-study scripts. Only two
things are replaced for the browser: multiprocessing (not available, runs
are executed serially) and SciPy's beta quantile in the Clopper-Pearson
interval (a pure-Python incomplete-beta inversion, checked against SciPy to
1e-10). The ellipse bootstrap uses the same statistic without the unused
Gaussian-radius term, which only affects speed.
"""
import json
import math
import random
import sys
import types

if "multiprocessing" not in sys.modules or not hasattr(sys.modules["multiprocessing"], "Pool"):
    _mp = types.ModuleType("multiprocessing")

    class _NoPool:
        def __init__(self, *a, **k):
            raise RuntimeError("multiprocessing is not available in the browser")
    _mp.Pool = _NoPool
    sys.modules["multiprocessing"] = _mp

import os  # noqa: E402

# Repository root: /repo inside the browser (Pyodide); the checkout folder when run natively.
ROOT = os.environ.get("SHYMA_ROOT") or ("/repo" if os.path.isfile("/repo/shyma/__init__.py")
                                        else os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for _p in (os.path.join(ROOT, "scripts"), ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run_casestudy as cs          # noqa: E402
import run_montecarlo as mc         # noqa: E402
import run_vsb30_validation as val  # noqa: E402
import run_nesc_checkcases as ne    # noqa: E402
import launch_optimizer as lo       # noqa: E402
from launch_optimizer import window_objects as _window_objects, diagnose as _diagnose  # noqa: E402,F401
from shyma import analysis, atmosphere, dispersion as D, earth, trajectory, vehicle  # noqa: E402

_PUB_FILE = os.path.join(ROOT, "web", "published.json")
if os.path.isfile(_PUB_FILE):          # embedded in the workbench by build_html.py
    PUB = json.load(open(_PUB_FILE))
else:                                  # native run: rebuild the published values from results/
    sys.path.insert(0, os.path.join(ROOT, "web"))
    import build_html as _bh           # noqa: E402
    PUB = _bh.published()
PAPER = {"window": dict(cs.WINDOW), "objective": dict(cs.OBJECTIVE), "steady": dict(cs.STEADY),
         "rn_le": cs.RN_LE, "beta_tv": cs.BETA_TV, "dsep": cs.DSEP, "coast_D": cs.COAST_D,
         "azimuth": cs.AZIMUTH, "nominal_el": dict(cs.NOMINAL), "payload": 400.0,
         "sigma": dict(D.SIG), "seed": mc.SEED, "n_seq": 1000, "phase": dict(cs.PHASE)}


# ---------------------------------------------------------------- utilities
def _clean(x):
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x


def _out(obj):
    return json.dumps(_clean(obj))


def _progress(kind, i, n, extra=None):
    try:
        from js import shymaProgress
        shymaProgress(json.dumps({"kind": kind, "i": i, "n": n, "extra": _clean(extra)}))
    except Exception:
        pass


def _betacf(a, b, x):
    tiny, qab, qap, qam = 1e-300, a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    return h


def betainc(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x)
    if x < (a + 1) / (a + b + 2):
        return math.exp(lb) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lb) * _betacf(b, a, 1 - x) / b


def beta_ppf(p, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if betainc(a, b, mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def clopper_pearson(k, n, conf=0.95):
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else beta_ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta_ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


D.clopper_pearson = clopper_pearson


def _set_windows(p):
    w = p["window"]
    cs.WINDOW = dict(mach=tuple(w["mach"]), alt_km=tuple(w["alt_km"]), q_kpa=tuple(w["q_kpa"]),
                     min_duration=float(w["min_duration"]))
    o = p["objective"]
    cs.OBJECTIVE = dict(mach=(float(o["mach_min"]), w["mach"][1]), alt_km=tuple(w["alt_km"]),
                        q_kpa=(float(o["q_min"]), w["q_kpa"][1]), min_duration=float(w["min_duration"]))
    cs.STEADY = dict(mach=tuple(w["mach"]), alt_km=tuple(w["alt_km"]), q_kpa=tuple(PAPER["steady"]["q_kpa"]),
                     min_duration=0.0)


def _restore():
    cs.WINDOW, cs.OBJECTIVE, cs.STEADY = (dict(PAPER["window"]), dict(PAPER["objective"]),
                                          dict(PAPER["steady"]))
    cs.RN_LE = PAPER["rn_le"]
    cs.PHASE.clear()
    cs.PHASE.update(PAPER["phase"])
    D.SIG.clear()
    D.SIG.update(PAPER["sigma"])


def window_polygon(w):
    """Mach-altitude outline of a window bounded by Mach, altitude and q (US76)."""
    left, right = [], []
    h0, h1 = w["alt_km"]
    n = 60
    for i in range(n + 1):
        h = h0 + (h1 - h0) * i / n
        p = atmosphere.atmosphere(h)[1]
        mq0 = math.sqrt(2 * w["q_kpa"][0] * 1e3 / (1.4 * p))
        mq1 = math.sqrt(2 * w["q_kpa"][1] * 1e3 / (1.4 * p))
        a, b = max(w["mach"][0], mq0), min(w["mach"][1], mq1)
        if a <= b:
            left.append([a, h])
            right.append([b, h])
    return left + right[::-1] + (left[:1] if left else [])




def _flight_heating(c, v, tr, s, dsep):
    """Peak leading-edge indicators after exposure, for runs without a test window."""
    t_min = v.t_burnout + dsep if c == "A" else None
    t_exp = t_min if c == "A" else 60.0
    if c == "D":
        t_exp = (s["t_window"] - 2.0) if s.get("t_window") is not None else v.t_burnout
    q, Q, Tw, qhw = analysis.heating_history(tr, cs.RN_LE, geometry="cylinder", t_expose=t_exp)
    ks = [k for k in range(len(tr.t)) if tr.t[k] >= t_exp and tr.alt[k] < 120e3] or [len(tr.t) - 1]
    kq = max(ks, key=lambda k: qhw[k])
    t0 = [analysis.T_from_h(analysis.h_air((tr.v[k] / tr.mach[k]) ** 2 / (1.4 * 287.05287)) + 0.5 * tr.v[k] ** 2)
          for k in ks if tr.mach[k] > 0]
    return {"t_expose": t_exp, "qdot_hw_peak_MW": qhw[kq] / 1e6, "t_peak": tr.t[kq],
            "mach_at_peak": tr.mach[kq], "alt_at_peak_km": tr.alt[kq] / 1e3,
            "T_w_peak_K": max(Tw[k] for k in ks), "T0_peak_K": max(t0) if t0 else None,
            "Q_to_impact_MJ": Q[-1] / 1e6}


# ---------------------------------------------------------------- API
def defaults():
    return _out({"paper": {**PAPER, "window": {"mach": PAPER["window"]["mach"], "alt_km": PAPER["window"]["alt_km"],
                                               "q_kpa": PAPER["window"]["q_kpa"],
                                               "min_duration": PAPER["window"]["min_duration"]},
                           "objective": {"mach_min": PAPER["objective"]["mach"][0],
                                         "q_min": PAPER["objective"]["q_kpa"][0]}},
                 "published": {"montecarlo_meta": PUB["montecarlo_meta"], "k_wind": PUB["k_wind"]}})


def _is_paper_nominal(p):
    c = p["concept"]
    ok = (abs(p["payload"] - 400.0) < 1e-9 and abs(p["el"] - PAPER["nominal_el"][c]) < 1e-9
          and abs(p["az"] - PAPER["azimuth"]) < 1e-9 and abs(p["beta"] - PAPER["beta_tv"]) < 1e-9
          and abs(p["rn_le"] - PAPER["rn_le"]) < 1e-12 and abs(p["dsep"] - PAPER["dsep"]) < 1e-9
          and (c != "D" or abs(p["coast"] - PAPER["coast_D"]) < 1e-9))
    w = p["window"]
    ok = ok and list(w["mach"]) == list(PAPER["window"]["mach"]) and list(w["alt_km"]) == list(PAPER["window"]["alt_km"]) \
        and list(w["q_kpa"]) == list(PAPER["window"]["q_kpa"]) and w["min_duration"] == PAPER["window"]["min_duration"] \
        and p["objective"]["mach_min"] == PAPER["objective"]["mach"][0] and p["objective"]["q_min"] == PAPER["objective"]["q_kpa"][0]
    return ok


def nominal(pjson):
    p = json.loads(pjson)
    c = p["concept"]
    try:
        _set_windows(p)
        cs.RN_LE = float(p["rn_le"])
        if p.get("phase", "default") != "default":
            cs.PHASE[c] = p["phase"]
        kw = {"coast": float(p["coast"])} if c == "D" else {}
        _progress("nominal", 0, 2)
        v, tr = cs.fly(c, float(p["payload"]), float(p["el"]), az=float(p["az"]),
                       beta_tv=float(p["beta"]), dt=0.01, rec=1, **kw)   # as run_casestudy.py (nominal)
        _progress("nominal", 1, 2)
        s = cs.summarize(c, v, tr, full=True, dsep=float(p["dsep"]))
        wins = _window_objects(c, v, float(p["dsep"]))
        s["diagnostics"] = {k: _diagnose(w, tr, c) for k, w in wins.items()}
        s["flight_heating"] = _flight_heating(c, v, tr, s, float(p["dsep"]))
        s["window_passes"] = _passes(wins["threshold"], tr)
        s["test_phase"] = cs.PHASE[c]
        poly = {"threshold": window_polygon(cs.WINDOW), "objective": window_polygon(cs.OBJECTIVE)}
        paper = _is_paper_nominal(p) and p.get("phase", "default") == "default"
        return _out({"summary": s, "published": PUB["nominal"][c] if paper else None,
                     "is_paper_case": paper, "polygons": poly})
    finally:
        _restore()


def _passes(w, tr):
    """Every crossing of the Mach-altitude-q box, with the reason it does or does not
    count as a test (flight phase, booster thrust, readiness time t_min)."""
    ta = tr.apogee()[0]
    out, start = [], None
    n = len(tr.t)

    def geo(k):
        a, m, q = tr.alt[k] / 1e3, tr.mach[k], tr.q[k] / 1e3
        return w.alt_km[0] <= a <= w.alt_km[1] and w.mach[0] <= m <= w.mach[1] and w.q_kpa[0] <= q <= w.q_kpa[1]
    for k in range(n):
        ok = geo(k)
        if ok and start is None:
            start = k
        if (not ok or k == n - 1) and start is not None:
            k1 = k - 1 if not ok else k
            ks = range(start, k1 + 1)
            ph = "ascent" if tr.t[k1] <= ta else "descent"
            reasons = []
            if w.phase not in ("any", ph):
                reasons.append("outside the test phase (%s)" % w.phase)
            if not w.powered and all(tr.thrust[i] > 0.0 for i in ks):
                reasons.append("booster thrusting")
            elif not w.powered and any(tr.thrust[i] > 0.0 for i in ks):
                reasons.append("booster thrusting during part of it")
            if w.t_min is not None and tr.t[k1] < w.t_min:
                reasons.append("before the vehicle is ready (t < %.1f s)" % w.t_min)
            elif w.t_min is not None and tr.t[start] < w.t_min:
                reasons.append("partly before the vehicle is ready (t < %.1f s)" % w.t_min)
            out.append({"t0": tr.t[start], "t1": tr.t[k1], "duration": tr.t[k1] - tr.t[start], "phase": ph,
                        "counted": not any(r.startswith("outside") or r == "booster thrusting" or r.startswith("before")
                                           for r in reasons),
                        "reasons": reasons})
            start = None
    return out




def optimize(pjson):
    p = json.loads(pjson)
    c = p["concept"]
    try:
        _set_windows(p)
        if p.get("phase", "default") != "default":
            cs.PHASE[c] = p["phase"]
        fixed = {"payload": float(p["payload"]), "az": float(p["az"]), "beta": float(p["beta"]),
                 "el": float(p["el"]), "coast": float(p["coast"]), "dsep": float(p["dsep"])}
        res = lo.search(c, fixed, p["free"], p["target"], p["objective_fn"], p.get("constraints", {}),
                        progress=lambda i, n, rec: _progress("opt", i, n, {"x": rec["x"], "dwell": rec[p["target"]],
                                                                            "feasible": rec["feasible"]}))
        bx = res["best_inputs"]
        po = PUB.get("optimizer", {})
        key = {("threshold", "max_dwell", 400.0): "max_dwell_400", ("objective", "max_dwell", 350.0): "objective_350"}.get(
            (p["target"], p["objective_fn"], fixed["payload"]))
        same_free = key and all(list(map(float, p["free"].get(k) or [])) == list(map(float, v))
                                for k, v in po.get("free", {}).get(c, {}).items()) and \
            all(not p["free"].get(k) for k in ("az", "dsep")) and (c == "D" or not p["free"].get("coast")) and \
            fixed["beta"] == PAPER["beta_tv"] and fixed["az"] == PAPER["azimuth"] and fixed["dsep"] == PAPER["dsep"] and \
            not p.get("constraints") and p.get("phase", "default") == "default"
        res["paper"] = po[key][c] if same_free else None
        res["final"] = json.loads(nominal(json.dumps(dict(p, payload=bx["payload"], el=bx["el"], az=bx["az"],
                                                          beta=bx["beta"], coast=bx["coast"], dsep=bx["dsep"]))))
        return _out(res)
    finally:
        _restore()


def _analyse(out, tau):
    n = len(out)
    lats = [o["impact_lat"] for o in out]
    lons = [o["impact_lon"] for o in out]
    e = D.ellipse(sum(lats) / n, sum(lons) / n, lats, lons)
    w = [o["window_s"] for o in out]
    met = sum(o["window_met"] for o in out)
    link = sum(1 for o in out if o["window_met"] and o.get("cla_margin_min") is not None
               and o.get("cla_el_min") is not None and o["cla_el_min"] >= getattr(mc, "ELEV_MASK_DEG", 3.0)
               and o["cla_margin_min"] - mc.ANT_LOSS_DB - mc.FADE_DB >= 0.0)
    obj = sum(o["objective_met"] for o in out)
    apo = [o["apogee_km"] for o in out]
    am = sum(apo) / n

    def smaj(sample):
        la = [a for a, _ in sample]
        lo = [b for _, b in sample]
        xs, ys = D._xy(sum(la) / len(la), sum(lo) / len(lo), la, lo)
        k = len(xs)
        mx, my = sum(xs) / k, sum(ys) / k
        sxx = sum((x - mx) ** 2 for x in xs) / (k - 1)
        syy = sum((y - my) ** 2 for y in ys) / (k - 1)
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (k - 1)
        t, det = sxx + syy, sxx * syy - sxy ** 2
        return math.sqrt(t / 2 + math.sqrt(max(t ** 2 / 4 - det, 0.0)))

    res = {"n": n, "p_window": met / n, "p_window_ci95": D.clopper_pearson(met, n),
           "p_window_link": link / n, "p_window_link_ci95": D.clopper_pearson(link, n),
           "p_objective": obj / n, "p_objective_ci95": D.clopper_pearson(obj, n),
           "window_p05": D.percentile(w, 0.05), "window_p50": D.percentile(w, 0.5),
           "window_p95": D.percentile(w, 0.95),
           "window_p05_ci95": D.bootstrap_ci(w, lambda x: D.percentile(x, 0.05)),
           "steady_p50": D.percentile([o["steady_s"] for o in out], 0.5),
           "apogee_mean": am,
           "apogee_std": math.sqrt(sum((a - am) ** 2 for a in apo) / (n - 1)) if n > 1 else None,
           "range_mean": sum(o["impact_range_km"] for o in out) / n,
           "ellipse": {k: v for k, v in e.items() if k not in ("x_km", "y_km")},
           "xy": [e["x_km"], e["y_km"]], "window_samples": w}
    res["ellipse"]["sigma_major_ci95"] = D.bootstrap_ci(list(zip(lats, lons)), smaj, n_boot=200)
    return res


def montecarlo(pjson):
    p = json.loads(pjson)
    c, n, seed = p["concept"], int(p["n"]), int(p["seed"])
    scale = float(p["scale"])
    try:
        if p.get("sigma"):
            D.SIG.update({k: float(v) for k, v in p["sigma"].items()})
        k = PUB["k_wind"] if p["wind"] else 0.0
        # every concept restarts the same aleatory stream (common random numbers, as run_montecarlo.run)
        rng = random.Random(seed)
        draws = [D.sample(rng, k, scale) for _ in range(n)]
        out = []
        for i, d in enumerate(draws):
            out.append(mc.one((c, d)))
            _progress("mc", i + 1, n, {"window_s": out[-1]["window_s"], "met": out[-1]["window_met"]})
        res = _analyse(out, cs.WINDOW["min_duration"])
        paper_like = (seed == PAPER["seed"] and scale == 1.0 and p["wind"] and D.SIG == PAPER["sigma"])
        cmp = None
        pub = PUB["montecarlo"][c]
        if paper_like:
            ref = pub["window_samples"][:n]
            diffs = [abs(a - b) for a, b in zip(res["window_samples"], ref)]
            cmp = {"n_compared": len(diffs), "max_abs_diff_s": max(diffs) if diffs else None,
                   "paper_prefix_p_window": sum(1 for x in ref if x >= cs.WINDOW["min_duration"]) / len(ref)}
        return _out({"result": res, "published": pub, "paper_like": paper_like, "comparison": cmp,
                     "runs": out})
    finally:
        _restore()


def vsb30():
    rows = val.load()
    recs = []
    for i, r in enumerate(rows):
        recs += val.held_out_set([r])
        _progress("vsb30", i + 1, len(rows), {"flight": r["flight"]})
    return _out({"records": recs, "summary": val.summary(recs), "published": PUB["vsb30"]})


def _ref_compare(tr, refs, times, key="alt_m", attr="alt"):
    errs, spreads = [], []
    for t in times:
        vals = [ne.interp(r["t"], r[key], t) for r in refs.values() if r["t"][0] <= t <= r["t"][-1]]
        if len(vals) < 2:
            continue
        vals.sort()
        med = vals[len(vals) // 2] if len(vals) % 2 else 0.5 * (vals[len(vals) // 2 - 1] + vals[len(vals) // 2])
        errs.append(ne.interp(tr.t, getattr(tr, attr), t) - med)
        spreads.append(vals[-1] - vals[0])
    return {"max_abs_err": max(abs(e) for e in errs), "max_spread": max(spreads),
            "final_err": errs[-1], "final_spread": spreads[-1], "n_ref": len(refs)}


def nesc():
    out = {}
    cases = (("case05", ne.case05), ("case09", lambda: ne.case_cannon(0.0, 1000.0)),
             ("case10", lambda: ne.case_cannon(1000.0, 0.0)))
    for i, (name, fn) in enumerate(cases):
        tr = fn()
        out[name] = {"alt": _ref_compare(tr, PUB["nesc"][name]["ref_profiles"], [0.5 * k for k in range(1, 61)]),
                     "profile": {"t": tr.t[::10], "alt_m": tr.alt[::10]},
                     "published": PUB["nesc"][name]["metrics"]}
        _progress("nesc", i + 1, 4, {"case": name})
    init = PUB["nesc"]["case17"]["init_row"]
    ne.load_case = lambda pattern: {"SIM 05": [init]}
    tr, info = ne.case17()
    tr.speed = [math.sqrt(a * a + b * b + c * c) for a, b, c in zip(tr.vn, tr.ve, tr.vd)]
    times = [20.0 + 2.0 * k for k in range(0, 58)] + [134.17]
    refs = PUB["nesc"]["case17"]["ref_profiles"]
    out["case17"] = {"alt": _ref_compare(tr, refs, times),
                     "speed": _ref_compare(tr, refs, times, key="v", attr="speed"),
                     "profile": {"t": tr.t[::5], "alt_m": tr.alt[::5]},
                     "published": PUB["nesc"]["case17"]["metrics"]}
    _progress("nesc", 4, 4, {"case": "case17"})
    return _out(out)


def checks():
    """Independent spot checks of the implemented formulas."""
    res = []

    def add(name, value, ref, tol, unit, source):
        res.append({"name": name, "value": value, "reference": ref, "tolerance": tol, "unit": unit,
                    "source": source, "pass": abs(value - ref) <= tol})
    T, p, rho, a = atmosphere.atmosphere(30.0)
    add("US76 temperature at 30 km", T, 226.509, 0.01, "K", "US Standard Atmosphere 1976, Table I")
    add("US76 pressure at 30 km", p, 1197.0, 1.0, "Pa", "US Standard Atmosphere 1976, Table I")
    add("US76 density at 30 km", rho, 0.018410, 2e-5, "kg/m^3", "US Standard Atmosphere 1976, Table I")
    T0, p0, r0, a0 = atmosphere.atmosphere(0.0)
    add("US76 sea-level speed of sound", a0, 340.294, 0.01, "m/s", "US Standard Atmosphere 1976")
    M = 7.0
    V = M * a
    add("Eq. (2): 1/2 rho V^2 = 1/2 gamma p M^2", 0.5 * rho * V * V, 0.5 * 1.4 * p * M * M, 1e-6, "Pa",
        "identity (V = M a, a^2 = gamma p / rho)")
    qsg = analysis.sutton_graves(rho, V, 0.002, "cylinder")
    add("Eq. (6): Sutton-Graves, 2 mm cylinder, Mach 7 at 30 km", qsg / 1e6,
        1.7415e-4 * math.sqrt(rho / 0.002) * V ** 3 / math.sqrt(2) / 1e6, 1e-9, "MW/m^2",
        "hand evaluation of k sqrt(rho/Rn) V^3 / sqrt(2)")
    tw, qh = analysis.wall_equilibrium(qsg, T, V, 0.85)
    add("Radiative balance eps sigma Tw^4 = q_hw", 0.85 * 5.670374419e-8 * tw ** 4 / 1e6, qh / 1e6, 1e-6,
        "MW/m^2", "Stefan-Boltzmann (CODATA 2018)")
    h1, h2 = analysis.h_air(1999.0), analysis.h_air(2001.0)
    cp = (h2 - h1) / 2.0
    add("Thermally perfect air: effective gamma at 2000 K", cp / (cp - 287.05287), 1.30, 0.01, "-",
        "paper text, Sec. IV-E (about 1.3)")
    earth.set_model("wgs84")
    g = earth.gravity(earth.WGS84_A, 0.0, 0.0)
    add("Gravity at the equator (point mass + J2, Eq. 3)", math.sqrt(sum(x * x for x in g)),
        3.986004418e14 / 6378137.0 ** 2 * (1 + 1.5 * 1.08263e-3), 1e-6, "m/s^2",
        "mu/a^2 (1 + 3/2 J2), WGS-84 constants")
    m = earth.link_margin_db(1.0e6)
    lam = 2.99792458e8 / 2.25e9
    ref = 3.0 - 20 * math.log10(4 * math.pi * 1.0e6 / lam) - 6.0 + 10.0 + 228.5991 - 10 * math.log10(2.0e6) - 11.0
    add("Eq. (S3): link margin at 1000 km slant range", m, ref, 1e-3, "dB",
        "SMAD link equation, hand evaluation (-10 log10 kB = 228.599 dB)")
    earth.reset_model()
    return _out(res)


# ------------------------------------------------------------------------------------------
# Geometry and stability derivatives (shyma.geometry, shyma.stability)
# ------------------------------------------------------------------------------------------
def geometry_defaults():
    """Estimated VSB-30 geometry (full stack) and its lift-off centre of gravity."""
    from shyma import geometry as G, vehicle as V
    g = G.vsb30_geometry()
    veh = V.vsb30(400.0)
    xcg = G.mass_properties(G.vsb30_mass_items(veh, g["stations"]), 0.0)[1]
    cfg = g["full"]
    return _out({"profile": cfg.profile, "d_ref": cfg.d_ref, "x_cg": xcg,
                 "fins": [{k: getattr(f, k) for k in ("name", "x_le", "c_root", "c_tip", "span", "x_tip_offset", "n")}
                          for f in cfg.fins],
                 "machs": [0.5, 0.9, 1.2, 1.5, 2, 3, 4, 5, 6, 7, 8]})


def stability(pjson):
    """Static and dynamic pitch derivatives of a user-defined body of revolution with fin sets."""
    from shyma import geometry as G, stability as St
    p = json.loads(pjson)
    fins = [G.FinSet(f.get("name", "F%d" % i), f["x_le"], f["c_root"], f["c_tip"], f["span"],
                     f["x_tip_offset"], int(f.get("n", 4))) for i, f in enumerate(p.get("fins", []))]
    cfg = G.Configuration(p["profile"], fins, float(p["d_ref"]))
    xcg = float(p["x_cg"])
    rows = []
    for M in p["machs"]:
        d = St.pitch_derivatives(cfg, float(M), xcg)
        rows.append({"mach": M, "cn_alpha": d["cn_alpha"], "x_cp": d["x_cp"], "margin_cal": d["margin_cal"],
                     "cm_alpha": d["cm_alpha"], "cmq": d["cmq"]})
    from shyma import drag as Dr
    t_c, r_le, ear, rem = float(p.get("t_c", 0.03)), float(p.get("r_le", 0.003)), float(p.get("exit_ratio", 0.72)), \
        float(p.get("re_per_m", 5e6))
    for f in fins:
        f.t_c = t_c
    for row in rows:
        c0 = Dr.components(cfg, float(row["mach"]), rem, False, r_le, ear)
        row.update({"cd_friction": c0["friction"], "cd_wave": c0["wave"], "cd_base": c0["base"], "cd_off": c0["total"],
                    "cd_on": Dr.components(cfg, float(row["mach"]), rem, True, r_le, ear)["total"]})
    comps = St.pitch_derivatives(cfg, float(p["machs"][0]), xcg)["components"]
    return _out({"rows": rows, "components": [{"name": n, "cn_alpha": c, "x_cp": x} for n, c, x in comps],
                 "length": cfg.length, "s_ref": cfg.s_ref, "volume": cfg.volume(), "wetted": cfg.wetted_area(),
                 "fins": [{"name": f.name, "area": f.area, "aspect": f.aspect,
                           "sweep_le_deg": math.degrees(f.sweep_le())} for f in fins]})
