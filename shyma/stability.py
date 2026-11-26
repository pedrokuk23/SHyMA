"""Static and dynamic pitch-stability derivatives of a finned body of revolution
from its geometry (shyma.geometry), by the classical component build-up used
at the conceptual stage (e.g. Walker and Wolowicz, NASA TM X-287, 1960, for the
X-15 at Mach 2-12):

* body: slender-body theory, segment by segment (dC_N/dalpha = 2 dS/S_ref per
  radian; centre of pressure of each conical segment as in Barrowman, 1967);
* fin panels: subsonic lifting-surface relation of Diederich (as used by
  Barrowman); supersonic linear theory, 4/beta with the tip correction
  1 - 1/(2 beta A) for supersonic leading edges and the slender-wing
  (Jones-Stewart) reduction for subsonic leading edges; at small angles of
  attack the hypersonic small-disturbance theory of Van Dyke, used in TM X-287,
  has the same initial slope; cubic blending across Mach 0.8-1.2;
* fin-body interference: slender-body factors K_W(B) and K_B(W) of Pitts,
  Nielsen and Kaattari (NACA Rep. 1307, 1957), K_W(B) + K_B(W) = (1 + tau)^2;
* pitch damping C_mq (per q d / 2V, referred to S_ref and d_ref): local
  angle-of-attack (strip) build-up, C_mq = -2 sum C_N_alpha,i ((x_i - x_cg)/d)^2,
  the approximation used for tail surfaces in TM X-287; the alpha-dot terms are
  neglected.

All derivatives are per radian, referred to S_ref = pi d_ref^2 / 4 and d_ref.
The results are linear (small angle of attack) and do not include shock
interference between components, transonic flow details, base-flow or jet
effects; they are intended for pre-Phase A screening.
"""
import bisect
import math

from . import geometry as G


def ellip_e(k):
    """Complete elliptic integral of the second kind E(k) (AGM method)."""
    if k >= 1.0:
        return 1.0
    a, b, s, p = 1.0, math.sqrt(1.0 - k * k), 0.5 * k * k, 0.5
    while abs(a - b) > 1e-15:
        a, b, c = 0.5 * (a + b), math.sqrt(a * b), 0.5 * (a - b)
        p *= 2.0
        s += p * c * c
    return math.pi / (2.0 * a) * (1.0 - s)


def panel_cna_subsonic(f, mach):
    """Normal-force slope of one panel (per rad, referred to the panel area), Diederich/Barrowman."""
    beta = math.sqrt(max(1.0 - mach * mach, 1e-6))
    A = f.aspect
    gc = math.cos(f.sweep_mid())
    return 2.0 * math.pi * A / (2.0 + math.sqrt(4.0 + (beta * A / gc) ** 2))


def panel_cna_supersonic(f, mach):
    """Normal-force slope of one panel (per rad, referred to the panel area), linear theory."""
    beta = math.sqrt(mach * mach - 1.0)
    A = f.aspect
    tip = 1.0 - 1.0 / (2.0 * beta * A) if beta * A >= 1.0 else 0.5 * beta * A
    c = 4.0 / beta * tip
    m = beta / math.tan(f.sweep_le()) if f.sweep_le() > 1e-9 else float("inf")
    if m < 1.0:                    # subsonic leading edge: slender-wing (Jones-Stewart) reduction
        c *= 0.5 * math.pi * m / ellip_e(math.sqrt(1.0 - m * m))
    return c


def _hermite(m, m0, m1, y0, y1):
    s = (m - m0) / (m1 - m0)
    return y0 + (y1 - y0) * s * s * (3.0 - 2.0 * s)


def panel_cna(f, mach, m_sub=0.8, m_sup=1.2):
    if mach <= m_sub:
        return panel_cna_subsonic(f, mach)
    if mach >= m_sup:
        return panel_cna_supersonic(f, mach)
    return _hermite(mach, m_sub, m_sup, panel_cna_subsonic(f, m_sub), panel_cna_supersonic(f, m_sup))


def panel_xcp(f, mach, m_sub=0.8, m_sup=1.2):
    """Axial centre of pressure of a panel: Barrowman (subsonic), area centroid (supersonic,
    uniform loading of linear theory), blended across Mach 0.8-1.2."""
    cr, ct, xt = f.c_root, f.c_tip, f.x_tip_offset
    x_sub = f.x_le + xt / 3.0 * (cr + 2 * ct) / (cr + ct) + (cr + ct - cr * ct / (cr + ct)) / 6.0
    x_sup = f.centroid_x()
    if mach <= m_sub:
        return x_sub
    if mach >= m_sup:
        return x_sup
    return _hermite(mach, m_sub, m_sup, x_sub, x_sup)


def k_wb(tau):
    """Slender-body lift ratio of the wing in presence of the body, K_W(B) (NACA Rep. 1307)."""
    if tau <= 0.0:
        return 1.0
    if tau >= 0.999:          # limit 2 (round-off of the closed form near tau = 1)
        return 2.0
    a = 1.0 / tau - tau
    num = (1 + tau ** 4) * (0.5 * math.atan(0.5 * a) + 0.25 * math.pi) - tau * tau * (a + 2.0 * math.atan(tau))
    return 2.0 / math.pi * num / (1.0 - tau) ** 2


def k_bw(tau):
    """Slender-body lift carried over onto the body, K_B(W) = (1 + tau)^2 - K_W(B)."""
    return (1.0 + tau) ** 2 - k_wb(tau)


def body_components(cfg):
    """Slender-body normal-force slope and centre of pressure of each body segment."""
    out = []
    sref_r2 = (cfg.d_ref / 2.0) ** 2
    x0, r0 = cfg.profile[0]
    if r0 > 0.0:                                   # blunt or truncated tip
        out.append(("body tip", 2.0 * r0 * r0 / sref_r2, x0))
    for xa, ra, xb, rb in cfg.segments():
        if abs(rb - ra) < 1e-12:
            continue
        cna = 2.0 * (rb * rb - ra * ra) / sref_r2
        L = xb - xa
        q = ra / rb
        xcp = xa + L / 3.0 * (1.0 + (1.0 - q) / (1.0 - q * q))
        out.append(("body %.2f-%.2f m" % (xa, xb), cna, xcp))
    return out


def fin_components(cfg, mach):
    out = []
    for f in cfg.fins:
        rb = cfg.radius(f.x_le + 0.5 * f.c_root)
        tau = rb / (rb + f.span)
        c1 = panel_cna(f, mach) * f.area / cfg.s_ref           # one panel, referred to S_ref
        cn = 0.5 * f.n * c1                                    # n equally spaced panels, pitch plane
        xcp = panel_xcp(f, mach)
        out.append(("fins %s" % f.name, cn * k_wb(tau), xcp))
        out.append(("body in presence of %s" % f.name, cn * k_bw(tau), xcp))
    return out


def static(cfg, mach):
    """Components [(name, C_N_alpha, x_cp)], total C_N_alpha (1/rad) and centre of pressure (m)."""
    comps = body_components(cfg) + fin_components(cfg, mach)
    cna = sum(c for _, c, _ in comps)
    xcp = sum(c * x for _, c, x in comps) / cna
    return {"components": comps, "cn_alpha": cna, "x_cp": xcp}


def pitch_derivatives(cfg, mach, x_cg):
    """C_N_alpha, C_m_alpha, static margin (calibers) and C_mq about the CG station x_cg."""
    s = static(cfg, mach)
    d = cfg.d_ref
    cma = -sum(c * (x - x_cg) / d for _, c, x in s["components"])
    cmq = -2.0 * sum(c * ((x - x_cg) / d) ** 2 for _, c, x in s["components"])
    return {"cn_alpha": s["cn_alpha"], "x_cp": s["x_cp"], "cm_alpha": cma,
            "margin_cal": (s["x_cp"] - x_cg) / d, "cmq": cmq, "components": s["components"]}


_TABLES = {}


class GeometricAero:
    """Mach- and time-dependent pitch derivatives for the pseudo-5-DOF ascent
    (shyma.trajectory5), from the geometry of each flight configuration and the
    mass items of the vehicle. configs[k] is the configuration after k stage
    drops. Component tables are precomputed on a Mach grid and interpolated."""

    def __init__(self, configs, mass_items, mach_grid=None):
        self.configs = configs
        self.items = mass_items
        self.d = configs[0].d_ref
        self.grid = mach_grid or [0.05 * i for i in range(0, 241)]
        key = (tuple(id(c) for c in configs), tuple(self.grid))
        if key not in _TABLES:
            _TABLES[key] = [[body_components(cfg) + fin_components(cfg, m) for m in self.grid] for cfg in configs]
        self.tab = _TABLES[key]
        self.lengths = [c.length for c in configs]

    def _config(self, veh, tc):
        k = sum(1 for s in veh.stages if tc >= s.t_drop)
        return min(k, len(self.configs) - 1)

    def components(self, k, mach):
        g = self.grid
        mach = min(max(mach, g[0]), g[-1])
        j = min(max(bisect.bisect_right(g, mach) - 1, 0), len(g) - 2)
        w = (mach - g[j]) / (g[j + 1] - g[j])
        a, b = self.tab[k][j], self.tab[k][j + 1]
        return [(ca[0], ca[1] + w * (cb[1] - ca[1]), ca[2] + w * (cb[2] - ca[2])) for ca, cb in zip(a, b)]

    def x_cg(self, t):
        return G.mass_properties(self.items, t)[1]

    def coeffs(self, veh, t, mach, tc=None):
        """(C_N_alpha, static margin in metres, C_mq) at time t and Mach number `mach`."""
        tc = t if tc is None else tc
        comps = self.components(self._config(veh, tc), mach)
        xcg = G.mass_properties(self.items, t, tc)[1]
        cna = sum(c for _, c, _ in comps)
        xcp = sum(c * x for _, c, x in comps) / cna
        cmq = -2.0 * sum(c * ((x - xcg) / self.d) ** 2 for _, c, x in comps)
        return cna, xcp - xcg, cmq

    def inertia(self, veh, t, tc=None):
        return G.mass_properties(self.items, t, tc)[2]

    def length(self, veh, t, tc=None):
        tc = t if tc is None else tc
        return self.lengths[self._config(veh, tc)]


_GEOM = {}


def vsb30_aero(veh, n_fins=3, t_c=0.03, fin_scale=1.0):
    """GeometricAero of a SHyMA VSB-30 vehicle (geometry cached per option set)."""
    key = (n_fins, t_c, fin_scale)
    if key not in _GEOM:
        _GEOM[key] = G.vsb30_geometry(n_fins, t_c, fin_scale)
    g = _GEOM[key]
    return GeometricAero([g["full"], g["upper"]], G.vsb30_mass_items(veh, g["stations"]))


# ---------------------------------------------------------------------------------------------
# Inverse check: derivatives from a damped pitch oscillation (cf. Greenberg, NACA TN 2340, 1951;
# Jaffe, JPL TR 32-757, 1965)
# ---------------------------------------------------------------------------------------------
def fit_damped_oscillation(t, y):
    """Least-squares fit y = exp(-lam t) (A cos(w t) + B sin(w t)) (Gauss-Newton from a
    zero-crossing / peak-ratio initial guess). Returns (lam, w, A, B, rms)."""
    zc = [t[i] - y[i] * (t[i + 1] - t[i]) / (y[i + 1] - y[i]) for i in range(len(y) - 1)
          if y[i] == 0.0 or y[i] * y[i + 1] < 0.0]
    w = math.pi / ((zc[-1] - zc[0]) / (len(zc) - 1)) if len(zc) > 1 else 1.0
    pk = [(t[i], abs(y[i])) for i in range(1, len(y) - 1) if abs(y[i]) >= abs(y[i - 1]) and abs(y[i]) >= abs(y[i + 1])]
    lam = math.log(pk[0][1] / pk[-1][1]) / (pk[-1][0] - pk[0][0]) if len(pk) > 1 and pk[-1][1] > 0 else 0.1
    A, B = y[0], 0.0
    p = [lam, w, A, B]
    for _ in range(50):
        JTJ = [[0.0] * 4 for _ in range(4)]
        JTr = [0.0] * 4
        for ti, yi in zip(t, y):
            e = math.exp(-p[0] * ti)
            c, s = math.cos(p[1] * ti), math.sin(p[1] * ti)
            f = e * (p[2] * c + p[3] * s)
            J = [-ti * f, e * ti * (-p[2] * s + p[3] * c), e * c, e * s]
            r = yi - f
            for a in range(4):
                JTr[a] += J[a] * r
                for b in range(4):
                    JTJ[a][b] += J[a] * J[b]
        dp = _solve4(JTJ, JTr)
        p = [p[i] + dp[i] for i in range(4)]
        if max(abs(dp[i]) / (abs(p[i]) + 1e-12) for i in range(4)) < 1e-12:
            break
    rms = math.sqrt(sum((yi - math.exp(-p[0] * ti) * (p[2] * math.cos(p[1] * ti) + p[3] * math.sin(p[1] * ti))) ** 2
                        for ti, yi in zip(t, y)) / len(t))
    return p[0], p[1], p[2], p[3], rms


def _solve4(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i]))
        M[i], M[piv] = M[piv], M[i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n + 1):
                M[r][c] -= f * M[i][c]
    x = [0.0] * n
    for i in reversed(range(n)):
        x[i] = (M[i][n] - sum(M[i][c] * x[c] for c in range(i + 1, n))) / M[i][i]
    return x


def derivatives_from_oscillation(lam, w, q, V, S, d, m, I, cn_alpha):
    """Invert the short-period model alpha'' + (dd + n) alpha' + (k + dd n) alpha = 0,
    n = q S C_N_alpha / (m V), k = -q S d C_m_alpha / I, dd = -q S d^2 C_mq / (2 V I),
    for C_m_alpha and C_mq given the fitted decay rate lam and damped frequency w."""
    n = q * S * cn_alpha / (m * V)
    dd = 2.0 * lam - n
    k = w * w + lam * lam - dd * n
    return {"cm_alpha": -k * I / (q * S * d), "cmq": -dd * 2.0 * V * I / (q * S * d * d)}
