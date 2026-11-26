"""Inverse design of the launch settings of a SHyMA concept.

Given the test vehicle (mass, ballistic coefficient), the test windows and
a goal, the free launch settings (elevation, azimuth, Concept D interstage
coast, Concept A burnout-to-test delay) are searched on a grid and refined
by a multi-start pattern search with step halving from the three best grid
points. Goals: maximise the dwell in the threshold or objective window,
maximise the q-band dwell, or minimise the impact range with the
window met; optional constraints on impact range and ascent dynamic
pressure. Infeasible trials are ranked by their distance to the window
(normalised margin of the closest eligible sample), so the search moves
towards the window and, when no setting meets it, reports the closest one
and the violated requirement.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs   # noqa: E402
from shyma import analysis   # noqa: E402

MINSTEP = {"el": 0.25, "coast": 0.5, "dsep": 0.25, "az": 1.0}


def window_objects(c, v, dsep):
    """Same windows as run_casestudy.windows(), returned as TestWindow objects."""
    t_min = v.t_burnout + dsep if c == "A" else None
    pw = c == "D"
    if pw:
        t_min = [st for st in v.stages if st.name == "S30"][0].t_ign
    ph = cs.PHASE[c]
    return {"threshold": analysis.TestWindow(phase=ph, t_min=t_min, powered=pw, **cs.WINDOW),
            "objective": analysis.TestWindow(phase=ph, t_min=t_min, powered=pw, **cs.OBJECTIVE)}


ASSUMED = {   # requirements the point-mass model does not evaluate (reported as "not evaluated")
    "all": [("Attitude at window entry", "assumed nose-forward at small angle of attack (point-mass test article)"),
            ("Inlet start and engine operation", "not modeled")],
    "A": [("Separation and shroud jettison at about 80 kPa", "modeled as a time delay")],
    "B": [("Separation and exo-atmospheric reorientation", "assumed successful")],
    "C": [("Exo-atmospheric reorientation of the stack", "assumed successful")],
    "D": [("Delayed S30 ignition, stack stability, shroud jettison", "assumed successful; stack fires along its axis")],
}


def add_assumed(out, concept):
    """Append the requirements that depend on an assumption, with status
    'not evaluated' (ok = None), so that the report separates what the model
    shows met, what it shows violated and what it does not evaluate."""
    for name, txt in ASSUMED["all"] + ASSUMED.get(concept, []):
        out["items"].append({"name": name, "ok": None, "status": "not evaluated", "text": txt})
    return out


def diagnose(w, tr, concept=None):
    """Which requirement is not met: eligibility (phase, readiness time, booster off),
    then the closest eligible sample to the Mach-altitude-q box. Each item has a
    status: 'met' or 'violated' by the model; with a concept, the items that
    rest on an assumption are added as 'not evaluated'."""
    ta = tr.apogee()[0]
    elig = []
    for k in range(len(tr.t)):
        if tr.thrust[k] > 0.0 and not w.powered:
            continue
        if w.phase == "ascent" and tr.t[k] > ta:
            continue
        if w.phase == "descent" and tr.t[k] < ta:
            continue
        if w.t_min is not None and tr.t[k] < w.t_min:
            continue
        elig.append(k)
    res = w.evaluate(tr)
    out = {"met": res["met"], "dwell_s": res["duration"], "min_dwell_s": w.min_duration,
           "phase": w.phase, "t_min": w.t_min, "powered": w.powered, "n_eligible": len(elig), "items": []}
    if not elig:
        out["items"].append({"name": "eligible phase", "ok": False,
                             "text": "no sample in the %s phase after t = %s s with the booster off"
                                     % (w.phase, "%.1f" % w.t_min if w.t_min is not None else "0"),
                             "status": "violated"})
        return add_assumed(out, concept) if concept else out

    def marg(k):
        a, m, q = tr.alt[k] / 1e3, tr.mach[k], tr.q[k] / 1e3
        return [(a - w.alt_km[0]), (w.alt_km[1] - a), 10 * (m - w.mach[0]), 10 * (w.mach[1] - m),
                0.1 * (q - w.q_kpa[0]), 0.1 * (w.q_kpa[1] - q)]
    kb = max(elig, key=lambda k: min(marg(k)))
    a, m, q = tr.alt[kb] / 1e3, tr.mach[kb], tr.q[kb] / 1e3
    em = [tr.mach[k] for k in elig]
    ea = [tr.alt[k] / 1e3 for k in elig]
    eq = [tr.q[k] / 1e3 for k in elig]
    out["closest"] = {"t": tr.t[kb], "mach": m, "alt_km": a, "q_kpa": q}
    out["range"] = {"mach": [min(em), max(em)], "alt_km": [min(ea), max(ea)], "q_kpa": [min(eq), max(eq)]}

    def item(name, val, lo, hi, unit, rng):
        ok = lo <= val <= hi
        if ok:
            txt = "%.2f%s inside [%g, %g]" % (val, unit, lo, hi)
        else:
            side = "below" if val < lo else "above"
            txt = "%.2f%s %s the %s bound %g%s (eligible flight spans %.2f to %.2f%s)" % (
                val, unit, side, "lower" if val < lo else "upper", lo if val < lo else hi, unit, rng[0], rng[1], unit)
        out["items"].append({"name": name, "ok": ok, "text": txt, "status": "met" if ok else "violated"})
    item("Mach", m, w.mach[0], w.mach[1], "", out["range"]["mach"])
    item("Altitude", a, w.alt_km[0], w.alt_km[1], " km", out["range"]["alt_km"])
    item("Dynamic pressure", q, w.q_kpa[0], w.q_kpa[1], " kPa", out["range"]["q_kpa"])
    out["items"].append({"name": "Dwell", "ok": res["met"],
                         "text": "%.2f s %s minimum %.1f s" % (res["duration"], ">=" if res["met"] else "<", w.min_duration),
                         "status": "met" if res["met"] else "violated"})
    return add_assumed(out, concept) if concept else out


def quick(c, payload, el, az, beta, coast, dsep, target, dt=0.05):
    """Fast evaluation used by the optimizer (no heating post-processing)."""
    kw = {"coast": coast} if c == "D" else {}
    v, tr = cs.fly(c, payload, el, az=az, beta_tv=beta, dt=dt, rec=2, **kw)
    wins = window_objects(c, v, dsep)
    r = {}
    for name, w in wins.items():
        e = w.evaluate(tr)
        r[name] = e["duration"]
        r[name + "_met"] = e["met"]
    st = analysis.TestWindow(phase=cs.PHASE[c], t_min=wins["threshold"].t_min, powered=wins["threshold"].powered,
                             **cs.STEADY).evaluate(tr) if r["threshold_met"] else {"duration": 0.0}
    r["steady"] = st["duration"]
    d = diagnose(wins[target], tr)
    if "closest" in d:
        w = wins[target]
        cl = d["closest"]
        r["closeness"] = min((cl["alt_km"] - w.alt_km[0]), (w.alt_km[1] - cl["alt_km"]),
                             10 * (cl["mach"] - w.mach[0]), 10 * (w.mach[1] - cl["mach"]),
                             0.1 * (cl["q_kpa"] - w.q_kpa[0]), 0.1 * (w.q_kpa[1] - cl["q_kpa"]))
    else:
        r["closeness"] = -100.0
    ta = tr.apogee()[0]
    r["impact_range_km"] = analysis.downrange_km(tr, cs.SITE)[-1]
    r["q_max_ascent_kpa"] = max(tr.q[k] for k in range(len(tr.t)) if tr.t[k] <= ta) / 1e3
    return r


def score(r, objective, target, cons):
    dwell = r[target]
    met = r[target + "_met"]
    viol = 0.0
    if cons.get("range_max") is not None:
        viol += max(0.0, r["impact_range_km"] - cons["range_max"]) / 50.0
    if cons.get("q_max") is not None:
        viol += max(0.0, r["q_max_ascent_kpa"] - cons["q_max"]) / 20.0
    feasible = met and viol == 0.0
    if objective == "max_dwell":
        base = dwell if dwell > 0 else min(r["closeness"], 0.0) - 1e-3
    elif objective == "max_steady":
        base = r["steady"] if met else (dwell if dwell > 0 else min(r["closeness"], 0.0) - 1e-3) - 100.0
    else:  # min_range: feasible points ranked by range
        base = (-r["impact_range_km"]) if met else (dwell if dwell > 0 else min(r["closeness"], 0.0) - 1e-3) - 1e4
    # feasible settings always rank above infeasible ones; the penalty only orders the infeasible ones
    return base - 1e3 * viol - (0.0 if feasible else 1e6), feasible



def search(c, fixed, free, target="threshold", objective="max_dwell", cons=None, progress=None, dt=0.05):
    """Grid plus multi-start pattern search. fixed: payload, az, beta, el, coast, dsep;
    free: {name: (min, max, grid step)} for el, az, coast (D) and dsep (A)."""
    cons = cons or {}
    free = {k: v for k, v in free.items() if v and (k != "coast" or c == "D") and (k != "dsep" or c == "A")}
    axes = {}
    for k, rg in free.items():
        lo, hi, st = float(rg[0]), float(rg[1]), max(float(rg[2]), 1e-6)
        n = int(round((hi - lo) / st)) + 1
        axes[k] = [lo + i * st for i in range(n) if lo + i * st <= hi + 1e-9]
    grid = [{}]
    for k, vals in axes.items():
        grid = [dict(g, **{k: x}) for g in grid for x in vals]
    n_eval = len(grid)
    evals = []

    def run(pt):
        q = dict(fixed, **pt)
        r = quick(c, q["payload"], q["el"], q["az"], q["beta"], q["coast"], q["dsep"], target, dt)
        sc, ok = score(r, objective, target, cons)
        rec = {"x": {k: q[k] for k in axes}, "score": sc, "feasible": ok, **r}
        evals.append(rec)
        if progress:
            progress(len(evals), n_eval, rec)
        return rec
    best = None
    for pt in grid:
        rec = run(pt)
        if best is None or rec["score"] > best["score"]:
            best = rec
    # local refinement: pattern search with step halving from the three best grid points
    lims = {k: (float(free[k][0]), float(free[k][1])) for k in axes}
    minstep = MINSTEP
    starts = sorted(evals, key=lambda e: -e["score"])[:3] if axes else []
    seen = {tuple(sorted(e["x"].items())) for e in evals}
    for st in starts:
        cur = st
        steps = {k: float(free[k][2]) / 2 for k in axes}
        it = 0
        while it < 30 and any(steps[k] >= minstep.get(k, 0.25) for k in axes):
            improved = False
            for k in axes:
                if steps[k] < minstep.get(k, 0.25):
                    continue
                for sgn in (1, -1):
                    x = dict(cur["x"])
                    x[k] = round(min(max(x[k] + sgn * steps[k], lims[k][0]), lims[k][1]), 6)
                    key = tuple(sorted(x.items()))
                    if key in seen:
                        continue
                    seen.add(key)
                    n_eval += 1
                    rec = run(x)
                    it += 1
                    if rec["score"] > cur["score"] + 1e-9:
                        cur, improved = rec, True
                        break
                if improved:
                    break
            if not improved:
                for k in axes:
                    steps[k] /= 2
        if cur["score"] > best["score"]:
            best = cur
    return {"best": best, "best_inputs": dict(fixed, **best["x"]), "feasible": best["feasible"],
            "any_feasible": any(e["feasible"] for e in evals), "evaluations": evals, "axes": list(axes),
            "n_eval": len(evals)}
