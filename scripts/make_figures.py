"""Figures for the manuscript (vector PDF, 90 mm single-column and 137 mm full text width of the
Elsevier preprint layout, so that they are included at their natural size).

Fonts are embedded as TrueType (Type 42), Times-metric (Liberation Serif),
with STIX math; minimum printed size is 7 pt at the natural width.
"""
import json
import sys
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, ROOT)
from shyma import dispersion as D  # noqa: E402
RES = os.environ.get("SHYMA_RESULTS", os.path.join(ROOT, "results"))
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

# colour-blind-safe categorical slots (Okabe-Ito-like)
C1, C2, C3, C4 = "#0072B2", "#D55E00", "#009E73", "#CC79A7"
INK, INK2, GRID = "#0b0b0b", "#4d4d4d", "#d9d8d4"
W1, W2 = 90 / 25.4, 137 / 25.4
# concept colours: four distinct Okabe-Ito hues (blue, vermillion, reddish purple, bluish green)
CCOL = {"A": C1, "B": C2, "C": C4, "D": C3}
# redundant line styles so the concepts stay distinguishable in grey-scale print
CLS = {"A": "-", "B": "--", "C": ":", "D": "-."}
WINCOL = "#7f7f7f"
TAU_VMAX = 15.0  # common colour scale of the dwell maps (s), saturated above  # neutral fill for test windows, so no concept colour is reused

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Liberation Serif", "Times New Roman", "Times"],
    "mathtext.fontset": "stix", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
    "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.unicode_minus": True,
})


def load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def save(fig, name):
    fig.savefig(os.path.join(FIG, name))
    plt.close(fig)


def fig_vsb30():
    """Residuals per flight (FT-10): SHyMA with the +-1 sigma Monte Carlo flight dispersion,
    operator nominal predictions, flights grouped by site, acceptance band of the apogee."""
    d = load("vsb30_validation.json")["flights"]
    d = sorted(d, key=lambda f: (f["site"] != "CLA", f["flight"]))
    x = list(range(len(d)))
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6))
    for k, (em, en, sk, lab, scale) in enumerate((
            ("apogee_err_model_pct", "apogee_err_nominal_pct", "mc_sigma_apogee_km", "Apogee error (%)", "pct"),
            ("range_err_model_km", "range_err_nominal_km", "mc_sigma_range_km", "Ground-range error (km)", "km"))):
        a = ax[k]
        if k == 0:
            a.axhspan(-5, 5, color=GRID, alpha=0.6, lw=0, label=r"$\pm$5% (RMS acceptance)")
        a.axhline(0, color=INK2, lw=0.6)
        first_held = True
        for i, f in enumerate(d):
            sig = f.get(sk, 0.0)
            if scale == "pct":
                sig = 100 * sig / f["apogee_actual"]
            held = not f["calibration"]
            lab_m = None
            if held and first_held:
                lab_m = r"SHyMA, held out ($\pm1\sigma$ flight dispersion)"
                first_held = False
            elif not held:
                lab_m = "SHyMA, calibration flight (V01)"
            a.errorbar([i - 0.13], [f[em]], yerr=[sig] if held else None, fmt="o" if held else "*",
                       ms=4.5 if held else 8, color=C1 if held else C3, mec="white", mew=0.5,
                       elinewidth=0.8, capsize=0, label=lab_m)
            a.plot([i + 0.13], [f[en]], "s", ms=4.5, mfc="none", mec=C2, mew=1.0,
                   label="Operator pre-flight nominal" if i == 0 else None)
        ncla = sum(f["site"] == "CLA" for f in d)
        a.axvline(ncla - 0.5, color=INK2, lw=0.6, ls=":")
        a.set_xticks(x)
        a.set_xticklabels([f["flight"] for f in d], rotation=90)
        a.set_ylabel(lab)
        a.grid(axis="x", visible=False)
        a.set_xlim(-0.6, len(d) - 0.4)
        yl = a.get_ylim()
        a.text(ncla / 2 - 0.5, yl[1], "CLA", ha="center", va="bottom", fontsize=7, color=INK2)
        a.text((ncla + len(d)) / 2 - 0.5, yl[1], "Esrange", ha="center", va="bottom", fontsize=7, color=INK2)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.0))
    ax[0].set_title("(a) Apogee", loc="left", pad=11)
    ax[1].set_title("(b) Ground range", loc="left", pad=11)
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig_vsb30_validation.pdf")


def fig_hyshot():
    d = load("hyshot_validation.json")
    p = d["profile"]
    fl = d["flight"]
    desc = [k for k in range(len(p["t"])) if p["t"][k] > 300 and 15 <= p["alt_km"][k] <= 60]
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6), sharey=True)
    ax[0].plot([p["mach"][k] for k in desc], [p["alt_km"][k] for k in desc], color=C1,
               label=r"SHyMA, $\beta = 5000$ kg/m$^2$")
    ax[1].plot([p["q_kpa"][k] for k in desc], [p["alt_km"][k] for k in desc], color=C1)
    sty = {"beta_1500": ("v", C3, r"$\beta = 1500$ kg/m$^2$ (35, 29 km)"),
           "beta_15000": ("^", C3, r"$\beta = 15000$ kg/m$^2$ (35, 29 km)"),
           "rho_x0.85": ("D", C4, r"density, pressure $\times 0.85$ (35, 29 km)")}
    for name, (mk, col, lab) in sty.items():
        s_ = d["sensitivity"][name]
        ax[0].plot([s_["mach_35"], s_["mach_29"]], [35, 29], ls="none", color=col, marker=mk, ms=4,
                   mfc="none", mew=1.0, label=lab)
        ax[1].plot([s_["q_35_kpa"], s_["q_29_kpa"]], [35, 29], ls="none", color=col, marker=mk, ms=4,
                   mfc="none", mew=1.0)
    ax[0].errorbar(fl["mach"], fl["alt_km"], xerr=0.2, fmt="s", ms=4.5, mfc="none", mec=C2, mew=1.2,
                   ecolor=C2, elinewidth=0.8, capsize=2, zorder=5,
                   label=r"HyShot 2 flight ($\pm$ acceptance tolerance)")
    ax[1].errorbar(fl["q_kpa"], fl["alt_km"], xerr=[0.1 * q for q in fl["q_kpa"]], fmt="s", ms=4.5,
                   mfc="none", mec=C2, mew=1.2, ecolor=C2, elinewidth=0.8, capsize=2, zorder=5)
    for a in ax:
        a.axhspan(29, 35, color=WINCOL, alpha=0.15, lw=0)
    allm = list(fl["mach"]) + [d["sensitivity"][n]["mach_%d" % h] for n in sty for h in (35, 29)]
    ax[0].set_xlim(min(allm) - 0.3, max(allm) + 0.3)
    allq = list(fl["q_kpa"]) + [d["sensitivity"][n]["q_%d_kpa" % h] for n in sty for h in (35, 29)]
    ax[1].set_xlim(0, 1.15 * max(allq))
    ax[0].set_xlabel("Mach number")
    ax[1].set_xlabel("Dynamic pressure (kPa)")
    ax[0].set_ylabel("Altitude (km)")
    ax[0].set_ylim(20, 45)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.02))
    ax[0].set_title("(a) Mach number", loc="left")
    ax[1].set_title("(b) Dynamic pressure", loc="left")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig_hyshot_validation.pdf")


def _window_band(ax, W, O):
    """True test window in the Mach-altitude plane: altitude, Mach and the
    dynamic-pressure bounds q = 0.7 p(h) M^2 (US76), shaded."""
    import sys
    sys.path.insert(0, ROOT)
    from shyma.atmosphere import atmosphere
    import numpy as np
    M = np.linspace(5.5, 8.5, 301)
    H = np.linspace(15, 45, 301)
    MM, HH = np.meshgrid(M, H)
    P = np.vectorize(lambda h: atmosphere(h)[1])(H)
    Q = 0.7 * P[:, None] * MM ** 2 / 1e3
    for win, col, alpha, lab in ((W, WINCOL, 0.20, "threshold window"), (O, WINCOL, 0.45, "objective window")):
        inside = ((MM >= win["mach"][0]) & (MM <= win["mach"][1]) & (HH >= win["alt_km"][0])
                  & (HH <= win["alt_km"][1]) & (Q >= win["q_kpa"][0]) & (Q <= win["q_kpa"][1]))
        ax.contourf(MM, HH, inside.astype(float), levels=[0.5, 1.5], colors=[col], alpha=alpha)
        ax.plot([], [], color=col, alpha=alpha + 0.2, lw=6, label=lab)


def fig_nominal():
    """(a) altitude vs downrange; (b) ascent-side tests A and D and (c) re-entry tests B and C in the
    Mach-altitude plane, with the in-window segment drawn thick in the concept colour (FT-02)."""
    import matplotlib.patheffects as pe
    d = load("casestudy.json")
    W, O = d["window"], d["objective"]
    fig = plt.figure(figsize=(W2, 4.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.5, wspace=0.08)
    ax0 = fig.add_subplot(gs[0, :])
    axb = fig.add_subplot(gs[1, 0])
    axc = fig.add_subplot(gs[1, 1], sharey=axb)
    N = d["nominal"]
    labs = {"A": "A: ascent test, %.0f°" % N["A"]["elevation"], "B": "B: re-entry, separated, %.0f°" % N["B"]["elevation"],
            "C": "C: re-entry, captive, %.0f°" % N["C"]["elevation"],
            "D": "D: delayed ignition, %.0f°, %.0f s coast" % (N["D"]["elevation"], d.get("coast_D", 24))}
    for c in "ABCD":
        n = N[c]
        p = n["profile"]
        col, ls = CCOL[c], CLS[c]
        ax0.plot(p["downrange_km"], p["alt_km"], color=col, ls=ls, lw=1.2, label=labs[c])
        ta = n["t_apogee"]
        up = c in "AD"
        ks_side = [k for k in range(len(p["t"])) if p["alt_km"][k] <= 60 and
                   ((p["t"][k] <= ta) if up else (p["t"][k] > ta))]
        a = axb if up else axc
        a.plot([p["mach"][k] for k in ks_side], [p["alt_km"][k] for k in ks_side], color=col, ls=ls, lw=1.2)
        if n.get("window_t") and n["window_t"][0] is not None:
            t0, t1 = n["window_t"]
            ks = [k for k in range(len(p["t"])) if t0 <= p["t"][k] <= t1]
            # A and B (drawn first) wide, D and C narrow on top, so that the coincident
            # in-window segments of A/D and B/C both stay visible
            wide = c in "AB"
            for aa, xk in ((ax0, "downrange_km"), (a, "mach")):
                aa.plot([p[xk][k] for k in ks], [p["alt_km"][k] for k in ks], color=col,
                        lw=5.0 if wide else 2.0, solid_capstyle="butt",
                        path_effects=[pe.Stroke(linewidth=6.6 if wide else 3.4, foreground="white"), pe.Normal()],
                        zorder=4 if wide else 5)
            k_lab = ks[0]
            off = {"A": (6, -12), "D": (-14, -4), "B": (6, 2), "C": (-14, -8)}[c]
            a.annotate(c, (p["mach"][k_lab], p["alt_km"][k_lab]), xytext=off,
                       textcoords="offset points", fontsize=8.5, color=col, fontweight="bold")
    for a, ttl in ((axb, "(b) Ascent side: A, D"), (axc, "(c) Re-entry: B, C")):
        _window_band(a, W, O)
        a.set_xlim(5.5, 8.5)
        a.set_ylim(15, 45)
        a.set_xlabel("Mach number")
        a.set_title(ttl, loc="left")
    axb.set_ylabel("Altitude (km)")
    plt.setp(axc.get_yticklabels(), visible=False)
    ax0.set_xlabel("Downrange from CLA (km)")
    ax0.set_ylabel("Altitude (km)")
    ax0.plot([], [], color=INK2, lw=4.0, label="time in test window (thick)")
    ax0.legend(loc="upper right", fontsize=7, ncol=1)
    h, l = axc.get_legend_handles_labels()
    axc.legend(h, l, loc="upper right", fontsize=7)
    ax0.set_title("(a) Trajectories", loc="left")
    save(fig, "fig_nominal_trajectories.pdf")


def fig_trade():
    import numpy as np
    cs_ = load("casestudy.json")
    d = cs_["grid"]
    els, pls = d["elevations"], d["payloads"]
    try:
        op = load("optimizer.json")["max_dwell_400"]
    except Exception:
        op = {}
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.5), sharey=True)
    titles = {"A": "(a) A: ascent", "B": "(b) B: separated", "C": "(c) C: captive"}
    for k, c in enumerate("ABC"):
        Z = np.array([[r["window_s"] for r in row] for row in d[c]])
        R = [[r["impact_range_km"] for r in row] for row in d[c]]
        im = ax[k].pcolormesh(els, pls, np.minimum(Z, TAU_VMAX), cmap="Blues", vmin=0, vmax=TAU_VMAX,
                              shading="nearest")
        cs = ax[k].contour(els, pls, R, levels=[200, 300, 400, 500, 600], colors=INK2, linewidths=0.6)
        ax[k].clabel(cs, cs.levels, fmt="%d", fontsize=6.5, inline_spacing=1)
        ax[k].contour(els, pls, Z, levels=[3.0], colors=C2, linewidths=1.6)
        if c in op:
            ax[k].plot([op[c]["x"]["el"]], [400], marker="x", ms=6, mew=1.5, color=INK)
        ax[k].set_xlabel("Elevation (deg)")
        ax[k].grid(False)
        ax[k].set_title(titles[c], loc="left")
    ax[0].set_ylabel("Test-vehicle mass (kg)")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, extend="max")
    cb.set_label("Time in threshold window (s)")
    h = [plt.Line2D([], [], color=C2, lw=1.6), plt.Line2D([], [], color=INK2, lw=0.6),
         plt.Line2D([], [], color=INK, marker="x", ls="none", ms=6, mew=1.5)]
    fig.legend(h, ["3 s requirement boundary", "impact range (km)", "longest dwell, 400 kg"],
               loc="upper center", ncol=3, bbox_to_anchor=(0.45, 0.0))
    save(fig, "fig_trade_space.pdf")


def fig_coverage():
    d = load("casestudy.json")["nominal"]
    fig, ax = plt.subplots(2, 1, figsize=(W1, 3.3))
    for k, c in enumerate("AB"):
        n = d[c]
        p = n["profile"]
        ax[k].plot(p["t"], p["el_CLA"], color=C1, label="CLA (Alcântara)")
        ax[k].plot(p["t"], p["el_CLBI"], color=C2, label="CLBI (Natal)")
        ax[k].axhline(3.0, color=INK2, lw=0.8, ls="--")
        t0, t1 = n["window_t"]
        ax[k].axvspan(t0, t1, color=C3, alpha=0.45, lw=0)
        ax[k].annotate("test window", xy=(0.5 * (t0 + t1), 60), xytext=(0.5 * (t0 + t1) - (120 if c == "B" else -40), 75),
                       fontsize=7, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
        ax[k].set_xlabel("Time from lift-off (s)")
        ax[k].set_ylim(-10, 90)
    ax[0].text(5, 5, "3° mask", color=INK2, fontsize=7)
    ax[0].set_ylabel("Elevation angle (deg)")
    ax[1].set_ylabel("Elevation angle (deg)")
    ax[0].legend(loc="upper right")
    ax[0].set_title("(a) Concept A", loc="left")
    ax[1].set_title("(b) Concept B (Concept C is similar)", loc="left")
    fig.tight_layout(h_pad=1.0)
    save(fig, "fig_coverage.pdf")


def fig_mc():
    d = load("montecarlo.json")["with_wind"]
    fig, axs = plt.subplots(2, 2, figsize=(W2, 5.2), sharex=True, sharey=True)
    ax = axs.ravel()
    titles = {"A": "(a) Concept A", "B": "(b) Concept B", "C": "(c) Concept C", "D": "(d) Concept D"}
    lim = 0.0
    for c in "ABCD":
        xs, ys = d[c]["xy"]
        mx, my = d[c]["ellipse"]["mean_xy"]
        lim = max(lim, max(abs(x - mx) for x in xs), max(abs(y - my) for y in ys))
    lim = 10 * math.ceil(lim / 10)
    for k, c in enumerate("ABCD"):
        e = d[c]["ellipse"]
        xs, ys = d[c]["xy"]
        mx, my = e["mean_xy"]
        ax[k].scatter([x - mx for x in xs], [y - my for y in ys], s=2.5, color=CCOL[c], alpha=0.5, lw=0)
        for nsig in (1, 3):
            ax[k].add_patch(Ellipse((0, 0), 2 * nsig * e["sigma_major_km"], 2 * nsig * e["sigma_minor_km"],
                                    angle=e["angle_from_east_deg"], fill=False, edgecolor=INK,
                                    lw=0.9, ls="-" if nsig == 3 else ":"))
        ax[k].add_patch(plt.Circle((0, 0), 50, fill=False, edgecolor=INK2, lw=1.0, ls="--"))
        ax[k].set_aspect("equal")
        ax[k].set_xlim(-lim, lim)
        ax[k].set_ylim(-lim, lim)
        ax[k].set_title(titles[c], loc="left")
    for a in axs[1]:
        a.set_xlabel("East (km)")
    for a in axs[:, 0]:
        a.set_ylabel("North (km)")
    h = [plt.Line2D([], [], color=INK, lw=0.9), plt.Line2D([], [], color=INK, lw=0.9, ls=":"),
         plt.Line2D([], [], color=INK2, lw=1.0, ls="--")]
    fig.legend(h, [r"$3\sigma$ ellipse", r"$1\sigma$ ellipse", "50 km VSB-30 design 3σ radius"],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.02))
    fig.tight_layout()
    save(fig, "fig_montecarlo.pdf")


def fig_mc_convergence():
    """(a) running P(tau_use >= 3 s) with the link, with its Clopper-Pearson 95 % band; (b) running
    5th-percentile window dwell with its order-statistic 95 % band (FT-05)."""
    from scipy.stats import binom
    d = load("montecarlo.json")["with_wind"]
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.5))
    pmin = []
    for c in "ABCD":
        u = [int(ch) for ch in d[c]["use_ok"]] if isinstance(d[c]["use_ok"], str) else [int(v) for v in d[c]["use_ok"]]
        w = d[c]["window_samples"]
        ns = list(range(50, len(w) + 1, 10))
        p = [sum(u[:n]) / n for n in ns]
        ci = [D.clopper_pearson(sum(u[:n]), n) for n in ns]
        pmin.append(min(x[0] for x in ci))
        ax[0].plot(ns, p, color=CCOL[c], ls=CLS[c], label="Concept " + c)
        ax[0].fill_between(ns, [x[0] for x in ci], [x[1] for x in ci], color=CCOL[c], alpha=0.12, lw=0)
        p05, l05, h05 = [], [], []
        for n in ns:
            srt = sorted(w[:n])
            p05.append(D.percentile(w[:n], 0.05))
            l05.append(srt[max(int(binom.ppf(0.025, n, 0.05)) - 1, 0)])
            h05.append(srt[min(int(binom.ppf(0.975, n, 0.05)), n - 1)])
        ax[1].plot(ns, p05, color=CCOL[c], ls=CLS[c], label="Concept " + c)
        ax[1].fill_between(ns, l05, h05, color=CCOL[c], alpha=0.12, lw=0)
    ax[0].set_xlabel("Number of runs")
    ax[1].set_xlabel("Number of runs")
    ax[0].set_ylabel(r"Running $P(\tau_{\rm use} \geq 3$ s$)$, link")
    ax[1].set_ylabel("Running 5th-pct. window dwell (s)")
    ax[0].set_ylim(max(0.0, 0.05 * math.floor(20 * min(pmin)) - 0.05), 1.005)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=4)
    ax[0].set_title("(a) Useful-time probability", loc="left")
    ax[1].set_title("(b) 5th-percentile dwell", loc="left")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig_mc_convergence.pdf")


def fig_robust():
    d = load("robust_selection.json")
    fig, ax = plt.subplots(1, 2, figsize=(W2, 2.6))
    MK = {"A": "o", "B": "s", "C": "^", "D": "D"}
    for c in "ABCD":
        r = d["concepts"][c]
        cf = r["confirm"]
        tags = [("nominal", ":", 1.0, "nominal")]
        if "primary" in cf and cf["primary"]["x"] != cf["nominal"]["x"]:
            tags.append(("primary", CLS[c], 1.4, "selected"))
        for tag, ls, lw, name in tags:
            cu = cf[tag]["curve"]
            el = cf[tag]["x"]["el"]
            sx = ("%.0f°" % el) + (("/%.0f s" % cf[tag]["x"]["coast"]) if c == "D" else "")
            if len(tags) == 1:
                ls, lw = CLS[c], 1.4
            ax[0].plot([x[0] for x in cu], [x[1] for x in cu], color=CCOL[c], ls=ls, lw=lw,
                       label="%s %s, %s" % (c, sx, name))
        xs = [q["range_mean_km"] for q in r["screen"]]
        ys = [q["p_req_ci95"][0] for q in r["screen"]]
        ax[1].scatter(xs, ys, s=12, marker=MK[c], facecolors="none", edgecolors=CCOL[c], lw=0.7,
                      alpha=0.8, label="Concept " + c + " candidates")
        nom = [q for q in r["screen"] if q["x"] == r["nominal_x"]][0]
        sel = [q for q in r["screen"] if q["x"] == r["primary_x"]][0]
        ax[1].scatter([nom["range_mean_km"]], [nom["p_req_ci95"][0]], s=40, marker=MK[c], color=CCOL[c],
                      edgecolors=INK, lw=0.6, zorder=5)
        ax[1].scatter([sel["range_mean_km"]], [sel["p_req_ci95"][0]], s=70, color=CCOL[c], marker="*",
                      edgecolors=INK, lw=0.4, zorder=6)
    ax[0].axvline(d["tau_req"], color=INK2, lw=0.8, ls="--")
    ax[0].set_xlabel(r"Required useful time $\tau_{\rm req}$ (s)")
    ax[0].set_ylabel(r"$P(\tau_{\rm use} \geq \tau_{\rm req})$, confirmation")
    ax[0].set_xlim(0, 18)
    ax[1].axhline(d["p_min"], color=INK2, lw=0.8, ls="--")
    ax[1].set_xlabel("Mean impact range (km)")
    ax[1].set_ylabel(r"Lower 95% bound of $P(\tau_{\rm use} \geq 3$ s$)$")
    ax[1].set_ylim(-0.02, 1.02)
    ax[0].set_title("(a) Useful test time", loc="left")
    ax[1].set_title("(b) Screened candidates", loc="left")
    h0, l0 = ax[0].get_legend_handles_labels()
    h1, l1 = ax[1].get_legend_handles_labels()
    h1 += [plt.Line2D([], [], marker="s", ls="none", color=INK2, mec=INK, ms=6),
           plt.Line2D([], [], marker="*", ls="none", color=INK2, mec=INK, ms=9)]
    l1 += ["nominal setting (filled)", "primary selection"]
    fig.legend(h0 + h1, l0 + l1, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.02), fontsize=7)
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig_robust_selection.pdf")


def fig_nesc():
    d = load("nesc_checkcases.json")
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.5))
    # (a) case 5
    c = d["case05"]
    first = True
    for k, r in c["ref_profiles"].items():
        ax[0].plot(r["t"], [a / 1e3 for a in r["alt_m"]], color=GRID, lw=2.4,
                   label="NESC reference simulations" if first else None)
        first = False
    p = c["profile"]
    st = max(1, len(p["t"]) // 20)
    ax[0].plot(p["t"][::st], [a / 1e3 for a in p["alt_m"]][::st], ls="none", marker="o", ms=3.5,
               mfc="none", mew=0.9, color=C3, label="SHyMA")
    ax[0].set_title("(a) Case 5: dropped sphere", loc="left")
    ax[0].legend(loc="lower left", fontsize=7)
    # (b) cases 9 and 10 (offset markers)
    first = True
    for case, col, mk in (("case09", C1, "o"), ("case10", C2, "^")):
        for k, r in d[case]["ref_profiles"].items():
            ax[1].plot(r["t"], [a / 1e3 for a in r["alt_m"]], color=GRID, lw=2.4,
                       label="NESC reference simulations" if first else None)
            first = False
        p = d[case]["profile"]
        st = max(1, len(p["t"]) // 16)
        off = 0 if case == "case09" else st // 2
        ax[1].plot(p["t"][off::st], [a / 1e3 for a in p["alt_m"]][off::st], ls="none", marker=mk, ms=3.5,
                   mfc="none", mew=0.9, color=col, label="SHyMA, case %d" % int(case[-2:]))
    ax[1].legend(loc="lower center", fontsize=7)
    ax[1].set_title("(b) Cases 9 and 10: cannonball", loc="left")
    # (c) case 17
    c = d["case17"]
    first = True
    for k, r in c["ref_profiles"].items():
        ax[2].plot(r["t"], [a / 1e3 for a in r["alt_m"]], color=GRID, lw=2.4,
                   label="NESC reference simulations" if first else None)
        first = False
    p = c["profile"]
    st = max(1, len(p["t"]) // 30)
    ax[2].plot(p["t"][::st], [a / 1e3 for a in p["alt_m"]][::st], ls="none", marker="o", ms=3.5,
               mfc="none", mew=0.9, color=C1, label="SHyMA (from t = 20 s)")
    ax[2].axvline(c["info"]["t_staging_s"], color=INK2, lw=0.8, ls="--")
    ax[2].text(c["info"]["t_staging_s"] + 3, 15, "2nd-stage\nignition", fontsize=7, color=INK2)
    ax[2].legend(loc="upper left", fontsize=7)
    ax[2].set_title("(c) Case 17: two-stage rocket", loc="left")
    for a in ax:
        a.set_xlabel("Time (s)")
    ax[0].set_ylabel("Altitude (km)")
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig_nesc.pdf")


def fig_wind():
    d = load("wind_calibration.json")
    ks = [float(k) for k in d["k_grid"]]
    tab = d["table"]
    g = [tab[str(k)]["r3_gauss"] for k in d["k_grid"]]
    e = [tab[str(k)]["r3_emp"] for k in d["k_grid"]]
    s_ = [3 * tab[str(k)]["range_sd"] for k in d["k_grid"]]
    fig, ax = plt.subplots(figsize=(W1, 2.5))
    ax.plot(ks, g, "-o", color=C1, ms=3, label=r"$3\sigma$ radius, Gaussian fit")
    ax.plot(ks, e, ":s", color=C1, ms=3, mfc="none", label=r"$3\sigma$ radius, empirical")
    ax.plot(ks, s_, "-^", color=C2, ms=3, label=r"$3\times$ downrange std. dev.")
    ax.axhline(50, color=C3, ls="--", lw=1.0, label="50 km VSB-30 design radius")
    lo, hi = d["k_wind_ci95"]
    ax.axvspan(lo, hi, color=C3, alpha=0.15, lw=0)
    ax.axvline(d["k_wind"], color=INK2, lw=0.8)
    ax.annotate(r"$k_w$ = %.2f m/s [%.2f, %.2f]" % (d["k_wind"], lo, hi), (d["k_wind"], 50),
                xytext=(10, -26), textcoords="offset points", fontsize=7, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))
    ax.set_xlabel(r"Residual-dispersion amplitude $k_w$ (m/s)")
    ax.set_ylabel("Impact dispersion (km)")
    ax.legend(loc="upper left", fontsize=7)
    save(fig, "fig_wind_calibration.pdf")


def fig_concept_d():
    import numpy as np
    g = load("casestudy.json")
    gd = g["grid"]["D"]
    els, cs_ = gd["elevations"], gd["coasts"]
    try:
        rsel = load("robust_selection.json")["concepts"]["D"]["primary_x"]
    except Exception:
        rsel = None
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.5), sharey=True)
    for k, p in enumerate((400, 350, 300)):
        rows = gd[str(p)]
        Z = np.array([[r["window_s"] for r in row] for row in rows])
        O = np.array([[1.0 if r["objective_met"] else 0.0 for r in row] for row in rows])
        im = ax[k].pcolormesh(els, cs_, np.minimum(Z, TAU_VMAX), cmap="Blues", vmin=0, vmax=TAU_VMAX,
                              shading="nearest")
        ax[k].contour(els, cs_, Z, levels=[3.0], colors=C2, linewidths=1.6)
        if O.max() > 0:
            ax[k].contourf(els, cs_, O, levels=[0.5, 1.5], colors="none", hatches=["////"])
            ax[k].contour(els, cs_, O, levels=[0.5], colors=INK, linewidths=0.8)
        if p == 400:
            ax[k].plot([g["nominal"]["D"]["elevation"]], [g["coast_D"]], marker="*", ms=10, color=C3,
                       mec="white", mew=0.6)
            if rsel:
                ax[k].plot([rsel["el"]], [rsel["coast"]], marker="D", ms=5.5, color=C4, mec="white", mew=0.6)
        ax[k].set_xlabel("Elevation (deg)")
        ax[k].set_title("(%s) %d kg" % ("abc"[k], p), loc="left")
        ax[k].grid(False)
    ax[0].set_ylabel("Interstage coast (s)")
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, extend="max")
    cb.set_label("Time in threshold window (s)")
    from matplotlib.patches import Patch
    h = [plt.Line2D([], [], color=C2, lw=1.6), Patch(facecolor="none", edgecolor=INK, hatch="////"),
         plt.Line2D([], [], color=C3, marker="*", ls="none", ms=9),
         plt.Line2D([], [], color=C4, marker="D", ls="none", ms=5.5)]
    fig.legend(h, ["3 s requirement", "objective window met", "nominal (400 kg)", "selected under uncertainty"],
               loc="upper center", ncol=2, bbox_to_anchor=(0.45, 0.0))
    save(fig, "fig_concept_d.pdf")


def fig_mass_ceiling():
    """Probability versus test-vehicle mass at the nominal settings (FT-19): observed grid proportions with
    Clopper-Pearson intervals and the pointwise one-sided 95 % lower band of the interpolated estimate."""
    d = load("mass_ceiling_mc.json")
    fig, axs = plt.subplots(2, 2, figsize=(W2, 4.4), sharey=True)
    for k, c in enumerate("ABCD"):
        a = axs.ravel()[k]
        r = d["concepts"][c]
        ms = [g_["mass_kg"] for g_ in r["grid"]]
        for key, col, lab, off in (("window", INK2, "window", -1.5), ("use", CCOL[c], "useful time, link", 1.5)):
            p = [g_[key]["p"] for g_ in r["grid"]]
            lo = [g_[key]["p"] - g_[key]["p_lo95"] for g_ in r["grid"]]
            hi = [g_[key]["p_hi95"] - g_[key]["p"] for g_ in r["grid"]]
            a.errorbar([m + off for m in ms], p, yerr=[lo, hi], fmt="o", ms=3, color=col,
                       mfc="white" if key == "window" else col, elinewidth=0.7, capsize=0,
                       label="observed, " + lab)
            a.plot(r[key]["curve_masses"], r[key]["curve_band"], color=col, lw=1.0,
                   ls="--" if key == "window" else "-", label="lower 95% band, " + lab)
        a.axhline(d["p_min"], color=INK2, lw=0.6, ls=":")
        a.axvline(400, color=INK2, lw=0.6, ls=":")
        a.set_title("(%s) Concept %s" % ("abcd"[k], c), loc="left")
        a.set_ylim(-0.02, 1.02)
    for a in axs[1]:
        a.set_xlabel("Test-vehicle mass (kg)")
    for a in axs[:, 0]:
        a.set_ylabel("Probability")
    h, l = axs[0, 1].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.02))
    fig.tight_layout()
    save(fig, "fig_mass_ceiling.pdf")


def _end_of_figures():
    pass


if __name__ == "__main__":
    FAILED = []
    for f in (fig_vsb30, fig_hyshot, fig_nominal, fig_trade, fig_coverage, fig_mc,
              fig_mc_convergence, fig_nesc, fig_wind, fig_concept_d, fig_robust, fig_mass_ceiling):
        try:
            f()
            print(f.__name__, "ok")
        except Exception as ex:  # report, continue, and fail at the end
            print(f.__name__, "FAILED", repr(ex))
            FAILED.append(f.__name__)
    if FAILED:
        sys.exit("figures failed: " + ", ".join(FAILED))
