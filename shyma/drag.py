"""Zero-lift drag of a finned body of revolution from its geometry (shyma.geometry),
by a conceptual-stage component build-up:

* skin friction: turbulent flat plate with compressibility by the Van Driest II
  method (adiabatic wall, recovery factor 0.89, viscosity power law T^0.76), on
  the body length and on the mean aerodynamic chord of the fins; subsonic form
  factors 1 + 60/(L/d)^3 + 0.0025 L/d (body) and 1 + 2 t/c (fins);
* body pressure (wave) drag at supersonic speed: shock-expansion method over the
  straight segments of the profile (Eggers and Savin, NACA Rep. 1249): the first
  segment by the exact Taylor-Maccoll cone solution, convex corners by
  Prandtl-Meyer expansion, concave corners by an oblique shock, both from the
  local state; a blunt tip as a hemispherical cap with the modified Newtonian
  pressure Cp_max cos^2;
* fin wave drag: two-dimensional shock-expansion on a diamond section in the plane
  normal to the leading edge (supersonic normal Mach number), c_d = Delta Cp_n
  cos^2(Lambda) t/c; leading-edge bluntness as a swept cylinder, (2/3) Cp_max
  cos^3(Lambda) on the frontal area;
* base drag: the conceptual-design correlation 0.12 + 0.13 M^2 (M < 1) and
  0.25/M (M > 1), referred to the base area, limited by the vacuum value
  2/(gamma M^2), and reduced by the factor (1 - A_exit/A_base) while a motor
  burns (as in the OpenRocket technical documentation, Niskanen 2009);
* transonic: the supersonic wave drags at Mach 1.5 are blended to zero at Mach 0.8
  with a cubic, and the form factors switch off between Mach 0.8 and 1.2.

All coefficients are referred to S_ref = pi d_ref^2 / 4. The model ignores
interference drag, protuberances, boattail separation and real-gas effects.
"""
import bisect
import math

from . import geometry as G

GAM = 1.4


# ----------------------------------------------------------------------------- gas dynamics
def pm_nu(M):
    g = GAM
    return math.sqrt((g + 1) / (g - 1)) * math.atan(math.sqrt((g - 1) / (g + 1) * (M * M - 1))) - \
        math.atan(math.sqrt(M * M - 1))


def pm_inverse(nu):
    lo, hi = 1.0, 60.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if pm_nu(mid) < nu:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def p_over_p0(M):
    return (1 + 0.5 * (GAM - 1) * M * M) ** (-GAM / (GAM - 1))


def oblique(M1, delta):
    """Weak oblique shock for deflection delta (rad): (beta, M2, p2/p1, p02/p01) or None if detached."""
    g = GAM
    mu = math.asin(1.0 / M1)
    best = None
    lo, hi = mu, 0.5 * math.pi
    # deflection as a function of beta has one maximum; find it, then solve on the weak branch
    def defl(b):
        return math.atan(2.0 / math.tan(b) * (M1 * M1 * math.sin(b) ** 2 - 1) /
                         (M1 * M1 * (g + math.cos(2 * b)) + 2))
    a, b_ = mu, 0.5 * math.pi
    for _ in range(80):
        m1, m2 = a + (b_ - a) / 3, b_ - (b_ - a) / 3
        if defl(m1) < defl(m2):
            a = m1
        else:
            b_ = m2
    bmax = 0.5 * (a + b_)
    if delta > defl(bmax):
        return None
    lo, hi = mu, bmax
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if defl(mid) < delta:
            lo = mid
        else:
            hi = mid
    beta = 0.5 * (lo + hi)
    mn1 = M1 * math.sin(beta)
    p21 = 1 + 2 * g / (g + 1) * (mn1 * mn1 - 1)
    mn2 = math.sqrt((1 + 0.5 * (g - 1) * mn1 * mn1) / (g * mn1 * mn1 - 0.5 * (g - 1)))
    M2 = mn2 / math.sin(beta - delta)
    p0r = ((g + 1) * mn1 * mn1 / ((g - 1) * mn1 * mn1 + 2)) ** (g / (g - 1)) * \
        ((g + 1) / (2 * g * mn1 * mn1 - (g - 1))) ** (1 / (g - 1))
    return beta, M2, p21, p0r


def cp_max(M):
    """Stagnation-point pressure coefficient behind a normal shock (modified Newtonian Cp_max)."""
    g = GAM
    if M <= 1.0:
        return (p_over_p0(M) ** -1 - 1) / (0.5 * g * M * M) if M > 0 else 1.0
    p02_p1 = ((g + 1) ** 2 * M * M / (4 * g * M * M - 2 * (g - 1))) ** (g / (g - 1)) * \
        (1 - g + 2 * g * M * M) / (g + 1)
    return (p02_p1 - 1) / (0.5 * g * M * M)


def taylor_maccoll(M, theta_c):
    """Exact inviscid flow over a sharp cone of half-angle theta_c (rad) at Mach M:
    returns (Cp_surface, M_surface, shock angle) or None if the shock is detached."""
    g = GAM
    k = 0.5 * (g - 1)

    def surface_angle(beta):
        r = oblique(M, 0.0) if False else None  # noqa: F841 (placeholder for clarity)
        mn1 = M * math.sin(beta)
        delta = math.atan(2.0 / math.tan(beta) * (mn1 * mn1 - 1) / (M * M * (g + math.cos(2 * beta)) + 2))
        mn2 = math.sqrt((1 + k * mn1 * mn1) / (g * mn1 * mn1 - k))
        M2 = mn2 / math.sin(beta - delta)
        V = (2.0 / ((g - 1) * M2 * M2) + 1) ** -0.5
        vr, vt = V * math.cos(beta - delta), -V * math.sin(beta - delta)
        th, h = beta, -1e-4

        def f(th_, vr_, vt_):
            a2 = k * (1 - vr_ * vr_ - vt_ * vt_)
            return vt_, (vt_ * vt_ * vr_ - a2 * (2 * vr_ + vt_ / math.tan(th_))) / (a2 - vt_ * vt_)
        while vt < 0.0 and th > 1e-4:
            k1 = f(th, vr, vt)
            k2 = f(th + 0.5 * h, vr + 0.5 * h * k1[0], vt + 0.5 * h * k1[1])
            k3 = f(th + 0.5 * h, vr + 0.5 * h * k2[0], vt + 0.5 * h * k2[1])
            k4 = f(th + h, vr + h * k3[0], vt + h * k3[1])
            nvr = vr + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            nvt = vt + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            if nvt >= 0.0:                       # interpolate the surface (v_theta = 0)
                w = vt / (vt - nvt)
                return th + w * h, vr + w * (nvr - vr), mn1
            th, vr, vt = th + h, nvr, nvt
        return th, vr, mn1
    mu = math.asin(1.0 / M)
    # the cone angle increases with the shock angle up to detachment
    lo, hi = mu + 1e-6, 0.5 * math.pi - 1e-6
    thc_hi = surface_angle(hi - 0.02)[0]
    if theta_c > thc_hi:
        return None
    hi = hi - 0.02
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if surface_angle(mid)[0] < theta_c:
            lo = mid
        else:
            hi = mid
    beta = 0.5 * (lo + hi)
    th, V, mn1 = surface_angle(beta)
    Ms = math.sqrt(2.0 / ((g - 1) * (1.0 / (V * V) - 1)))
    p0r = ((g + 1) * mn1 * mn1 / ((g - 1) * mn1 * mn1 + 2)) ** (g / (g - 1)) * \
        ((g + 1) / (2 * g * mn1 * mn1 - (g - 1))) ** (1 / (g - 1))
    p_pinf = p0r * p_over_p0(Ms) / p_over_p0(M)
    return (p_pinf - 1) / (0.5 * g * M * M), Ms, beta


# ----------------------------------------------------------------------------- friction
def cf_incompressible(re):
    """Karman-Schoenherr mean turbulent skin friction (Newton iteration)."""
    cf = 0.074 / re ** 0.2
    for _ in range(30):
        f = 0.242 / math.sqrt(cf) - math.log10(re * cf)
        df = -0.121 * cf ** -1.5 - 1.0 / (cf * math.log(10))
        cf -= f / df
    return cf


def cf_van_driest(re, M, r=0.89, n=0.76):
    """Van Driest II mean turbulent skin friction of an adiabatic flat plate."""
    if M < 0.1:
        return cf_incompressible(re)
    m = 0.5 * (GAM - 1) * M * M
    tw = 1 + r * m                       # adiabatic wall temperature / edge temperature
    A = math.sqrt(r * m / tw)
    B = (1 + r * m) / tw - 1
    den = math.sqrt(B * B + 4 * A * A)
    c1 = (2 * A * A - B) / den
    c2 = B / den
    fac = (math.asin(c1) + math.asin(c2)) / A
    cf = cf_incompressible(re) / tw ** 0.5
    for _ in range(60):
        lhs = 0.242 * fac / math.sqrt(cf * tw)
        rhs = math.log10(re * cf) - 0.5 * (1 + 2 * n) * math.log10(tw)
        g_ = lhs - rhs
        dg = -0.5 * lhs / cf - 1.0 / (cf * math.log(10))
        step = g_ / dg
        cf = max(cf - step, 1e-6)
        if abs(step) < 1e-12:
            break
    return cf


# ----------------------------------------------------------------------------- components
def _blend(M, m0=0.8, m1=1.2):
    s = min(max((M - m0) / (m1 - m0), 0.0), 1.0)
    return s * s * (3 - 2 * s)


def body_wave(cfg, M):
    """Supersonic pressure drag of the body forebody (shock-expansion), referred to S_ref."""
    if M <= 1.0:
        return 0.0
    q_inf = 0.5 * GAM * M * M
    cd = 0.0
    segs = cfg.segments()
    x0, r0 = cfg.profile[0]
    if r0 > 0.0:                                         # hemispherical cap, modified Newtonian
        cd += 0.5 * cp_max(M) * r0 * r0 / (cfg.d_ref / 2) ** 2
    Ml, pl, th_prev = None, None, None                   # local Mach, p/p_inf, previous angle
    for xa, ra, xb, rb in segs:
        th = math.atan2(rb - ra, xb - xa)
        if Ml is None:
            if th <= 1e-9:
                Ml, pl = M, 1.0
            else:
                tm = taylor_maccoll(M, th)
                if tm is None:                           # detached: Newtonian on the segment
                    cp = cp_max(M) * math.sin(th) ** 2
                    pl = 1 + cp * q_inf
                    Ml = max(1.05, M * math.cos(th))
                else:
                    cp, Ml, _ = tm
                    pl = 1 + cp * q_inf
        else:
            dth = th - th_prev
            if dth < -1e-9:                              # convex corner: expansion
                M2 = pm_inverse(pm_nu(Ml) - dth)
                pl *= p_over_p0(M2) / p_over_p0(Ml)
                Ml = M2
            elif dth > 1e-9:                             # concave corner: oblique shock
                ob = oblique(Ml, dth)
                if ob is None:
                    pl *= 1 + cp_max(Ml) * math.sin(dth) ** 2 * 0.5 * GAM * Ml * Ml
                    Ml = max(1.05, Ml * math.cos(dth))
                else:
                    pl *= ob[2]
                    Ml = ob[1]
        th_prev = th
        cp = (pl - 1) / q_inf
        cd += cp * (rb * rb - ra * ra) / (cfg.d_ref / 2) ** 2
    return cd


def fin_wave(f, M):
    """Supersonic wave drag of one fin panel (diamond section, sweep theory), referred to the panel area."""
    lam = f.sweep_le()
    Mn = M * math.cos(lam)
    if Mn <= 1.05 or f.t_c <= 0.0:
        return 0.0
    eps = math.atan(f.t_c / math.cos(lam))
    ob = oblique(Mn, eps)
    if ob is None:
        return 0.0
    _, M2, p21, _ = ob
    M3 = pm_inverse(pm_nu(M2) + 2 * eps)
    p31 = p21 * p_over_p0(M3) / p_over_p0(M2)
    dcp = (p21 - p31) / (0.5 * GAM * Mn * Mn)
    return dcp * math.cos(lam) ** 2 * f.t_c


def base_cd(M, powered=False, exit_area_ratio=0.0):
    """Base drag referred to the base area (power-off correlation, jet reduction)."""
    cd = 0.12 + 0.13 * M * M if M < 1.0 else 0.25 / M
    if M > 1.0:
        cd = min(cd, 2.0 / (GAM * M * M))
    return cd * (1.0 - exit_area_ratio) if powered else cd


class DragModel:
    """C_D(M, Re, powered) of the attached stack of each flight configuration."""

    def __init__(self, configs, r_le=0.003, exit_area_ratio=0.72, mach_grid=None):
        self.configs = configs
        self.r_le = r_le
        self.ear = exit_area_ratio
        self.grid = mach_grid or ([0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2]
                                  + [1.3 + 0.1 * i for i in range(8)] + [2.25 + 0.25 * i for i in range(32)])
        self.logre = [5.0 + 0.25 * i for i in range(21)]          # log10 Re on the body length
        self.wave = [[self._wave(cfg, M) for M in self.grid] for cfg in configs]
        self.fric = [[[self._friction(cfg, M, 10 ** lr) for lr in self.logre] for M in self.grid] for cfg in configs]

    def _wave(self, cfg, M):
        # below Mach 1.5 the shock-expansion method is unreliable (near-sonic local flow): the
        # supersonic wave drag at Mach 1.5 is blended to zero at Mach 0.8
        b = _blend(M, 0.8, 1.5)
        if M <= 0.8:
            return 0.0
        Ms = max(M, 1.5)
        w = body_wave(cfg, Ms)
        for f in cfg.fins:
            w += f.n * f.area / cfg.s_ref * fin_wave(f, Ms)
            lam = f.sweep_le()
            w += f.n * 2 * self.r_le * f.span / math.cos(lam) * math.cos(lam) / cfg.s_ref * \
                (2.0 / 3.0) * cp_max(Ms) * math.cos(lam) ** 3
        return b * w

    def _friction(self, cfg, M, re_l):
        L, d = cfg.length, 2 * max(r for _, r in cfg.profile)
        sub = 1.0 - _blend(M)
        ff_b = 1 + sub * (60.0 / (L / d) ** 3 + 0.0025 * L / d)
        cd = cf_van_driest(re_l, M) * cfg.wetted_area() / cfg.s_ref * ff_b
        for f in cfg.fins:
            mac = f.mac()[0]
            cd += cf_van_driest(re_l * mac / L, M) * 2 * f.n * f.area / cfg.s_ref * (1 + 2 * f.t_c)
        return cd

    def components(self, k, M, re_per_m, powered):
        cfg = self.configs[k]
        g = self.grid
        M = min(max(M, 0.0), g[-1])
        j = min(max(bisect.bisect_right(g, M) - 1, 0), len(g) - 2)
        w = (M - g[j]) / (g[j + 1] - g[j])
        wave = self.wave[k][j] + w * (self.wave[k][j + 1] - self.wave[k][j])
        lr = math.log10(max(re_per_m * cfg.length, 1e3)) if re_per_m else 7.5
        lr = min(max(lr, self.logre[0]), self.logre[-1])
        i = min(int((lr - self.logre[0]) / 0.25), len(self.logre) - 2)
        u = (lr - self.logre[i]) / 0.25
        fa = self.fric[k][j][i] + u * (self.fric[k][j][i + 1] - self.fric[k][j][i])
        fb = self.fric[k][j + 1][i] + u * (self.fric[k][j + 1][i + 1] - self.fric[k][j + 1][i])
        fric = fa + w * (fb - fa)
        base = base_cd(M, powered, self.ear) * (cfg.base_radius() / (cfg.d_ref / 2)) ** 2
        return {"friction": fric, "wave": wave, "base": base, "total": fric + wave + base}

    def cd(self, k, M, re_per_m=None, powered=False):
        return self.components(k, M, re_per_m, powered)["total"]


_MODELS = {}


def vsb30_drag(n_fins=3, t_c=0.03, fin_scale=1.0, r_le=0.003, exit_area_ratio=0.72):
    key = (n_fins, t_c, fin_scale, r_le, exit_area_ratio)
    if key not in _MODELS:
        g = G.vsb30_geometry(n_fins, t_c, fin_scale)
        _MODELS[key] = DragModel([g["full"], g["upper"]], r_le, exit_area_ratio)
    return _MODELS[key]


def attach(veh, model):
    """Use the geometric drag model for the attached stack of a vehicle.Vehicle
    (the payload phase keeps its own drag); cd_scale still multiplies it."""
    def cd_model(t, mach, tc, re_per_m):
        k = min(sum(1 for s in veh.stages if tc >= s.t_drop), len(model.configs) - 1)
        powered = any(s.burning(tc) and tc < s.t_drop for s in veh.stages)
        return model.cd(k, mach, re_per_m, powered)
    veh.cd_model = cd_model
    return veh


def components(cfg, M, re_per_m=5e6, powered=False, r_le=0.003, exit_area_ratio=0.72):
    """Drag components of one configuration at one Mach number, computed directly (no table)."""
    m = DragModel.__new__(DragModel)
    m.r_le, m.ear = r_le, exit_area_ratio
    wave = m._wave(cfg, M)
    fric = m._friction(cfg, M, re_per_m * cfg.length)
    base = base_cd(M, powered, exit_area_ratio) * (cfg.base_radius() / (cfg.d_ref / 2)) ** 2
    return {"friction": fric, "wave": wave, "base": base, "total": fric + wave + base}
