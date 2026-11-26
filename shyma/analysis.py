"""Mission-analysis post-processing for SHyMA: hypersonic test window,
stagnation-point aeroheating, tracking/telemetry coverage and impact point."""
import math

from . import earth
from .atmosphere import R_AIR

SUTTON_GRAVES_K = 1.7415e-4   # Earth air, SI units (W/m^2, kg/m^3, m, m/s)
SIGMA = 5.670374e-8


CYL_FACTOR = 1.0 / math.sqrt(2.0)   # 2-D (cylindrical) stagnation line vs sphere

# Thermally perfect air: translation + rotation (3.5 R T) plus harmonic
# vibration of N2 and O2. Dissociation is neglected (< 0.5 % of O2 at the
# stagnation conditions of the test windows).
_VIB = ((0.79, 3395.0), (0.21, 2239.0))   # (mole fraction, theta_v [K])


def h_air(T):
    """Specific enthalpy of thermally perfect air (J/kg), h(0) = 0."""
    hv = sum(x * th / math.expm1(th / T) for x, th in _VIB) if T > 1.0 else 0.0
    return R_AIR * (3.5 * T + hv)


def T_from_h(h, lo=50.0, hi=6000.0):
    """Inverse of h_air by bisection."""
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if h_air(m) > h:
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def sutton_graves(rho, v, rn, geometry="sphere"):
    """Cold-wall stagnation-point convective heat flux (W/m^2) for a nose
    (sphere) or leading-edge (unswept cylinder) radius rn (m)."""
    q = SUTTON_GRAVES_K * math.sqrt(rho / rn) * v ** 3
    return q * CYL_FACTOR if geometry == "cylinder" else q


def wall_equilibrium(q_cw, T_inf, v, emissivity=0.85):
    """Radiative-equilibrium wall temperature with the hot-wall correction
    q = q_cw (1 - h_w/h0), thermally perfect air. Returns (Tw, q_hw).
    Tw never exceeds the total temperature T0."""
    h0 = h_air(T_inf) + 0.5 * v * v
    T0 = T_from_h(h0)
    if q_cw <= 0.0:
        return T_inf, 0.0
    lo, hi = 0.0, T0
    for _ in range(60):
        Tw = 0.5 * (lo + hi)
        if emissivity * SIGMA * Tw ** 4 > q_cw * (1.0 - h_air(Tw) / h0):
            hi = Tw
        else:
            lo = Tw
    Tw = 0.5 * (lo + hi)
    return Tw, q_cw * (1.0 - h_air(Tw) / h0)


def wall_equilibrium_cpg(q_cw, T_inf, v, emissivity=0.85, cp=1004.5):
    """Same as wall_equilibrium for a calorically perfect gas (gamma = 1.4);
    kept only to quantify the effect of the thermally perfect model."""
    h0 = cp * T_inf + 0.5 * v * v
    T0 = h0 / cp
    if q_cw <= 0.0:
        return T_inf, 0.0
    lo, hi = 0.0, T0
    for _ in range(60):
        Tw = 0.5 * (lo + hi)
        if emissivity * SIGMA * Tw ** 4 > q_cw * (1.0 - Tw / T0):
            hi = Tw
        else:
            lo = Tw
    Tw = 0.5 * (lo + hi)
    return Tw, q_cw * (1.0 - Tw / T0)


def heating_history(tr, rn, emissivity=0.85, geometry="sphere", t_expose=None, atm=None):
    """Heat flux, heat load and wall temperature histories.

    q_cw : cold-wall flux; q_hw : hot-wall flux at radiative equilibrium;
    Tw   : radiative-equilibrium wall temperature (upper bound for a short,
           conduction-dominated exposure); Q : hot-wall heat load integrated
    from t_expose (shroud jettison / separation), zero before.
    """
    from .atmosphere import atmosphere as us76
    f = atm if atm is not None else us76
    q, qhw, Tw = [], [], []
    for r, v, a in zip(tr.rho, tr.v, tr.alt):
        qc = sutton_graves(r, v, rn, geometry) if a < 120e3 else 0.0
        T_inf = f(max(a, 0.0) / 1000.0)[0]
        tw, qh = wall_equilibrium(qc, T_inf, v, emissivity)
        q.append(qc)
        qhw.append(qh)
        Tw.append(tw)
    Q = [0.0]
    for k in range(1, len(q)):
        t0, t1 = tr.t[k - 1], tr.t[k]
        dq = 0.5 * (qhw[k] + qhw[k - 1]) * (t1 - t0)
        if t_expose is not None:
            if t1 <= t_expose:
                dq = 0.0
            elif t0 < t_expose and t1 > t0:   # count only the part after exposure
                dq *= (t1 - t_expose) / (t1 - t0)
        Q.append(Q[-1] + dq)
    return q, Q, Tw, qhw


class TestWindow:
    """Hypersonic test-window requirement (a priori stakeholder requirement).

    A sample is inside the window when Mach number, altitude and dynamic
    pressure are within bounds, the vehicle is unpowered, it is in the
    requested flight phase and t >= t_min (e.g. separation, shroud jettison
    and settling completed). Entry and exit times are located by linear
    interpolation of the constraint margins between samples, so the dwell
    time is not quantised by the output step.
    """

    def __init__(self, mach=(6.0, 8.0), alt_km=(23.0, 35.0), q_kpa=(20.0, 90.0),
                 min_duration=3.0, phase="descent", name="W", t_min=None, powered=False):
        self.mach = mach
        self.alt_km = alt_km
        self.q_kpa = q_kpa
        self.min_duration = min_duration
        self.phase = phase  # 'ascent', 'descent' or 'any'
        self.name = name
        self.t_min = t_min
        self.powered = powered   # True: boosted test, booster thrust allowed (F_b >= 0)

    def margin(self, tr, k, t_apogee):
        """Signed margin (>= 0 inside). Continuous in the flight variables."""
        if tr.thrust[k] > 0.0 and not self.powered:
            return -1.0
        a = tr.alt[k] / 1000.0
        q = tr.q[k] / 1000.0
        m = tr.mach[k]
        g = [(a - self.alt_km[0]), (self.alt_km[1] - a),          # km
             10.0 * (m - self.mach[0]), 10.0 * (self.mach[1] - m),  # 0.1 Mach
             0.1 * (q - self.q_kpa[0]), 0.1 * (self.q_kpa[1] - q)]  # 10 kPa
        if self.phase == "ascent":
            g.append(t_apogee - tr.t[k])
        elif self.phase == "descent":
            g.append(tr.t[k] - t_apogee)
        if self.t_min is not None:
            g.append(tr.t[k] - self.t_min)
        return min(g)

    def inside(self, tr, k, t_apogee):
        return self.margin(tr, k, t_apogee) >= 0.0

    def evaluate(self, tr):
        ta, _, _ = tr.apogee()
        g = [self.margin(tr, k, ta) for k in range(len(tr.t))]
        best, best_start, best_end = 0.0, None, None
        start = tr.t[0] if g[0] >= 0 else None
        for k in range(1, len(tr.t)):
            g0, g1 = g[k - 1], g[k]
            if g0 < 0 <= g1:      # entry
                w = -g0 / (g1 - g0) if (g1 > g0 and g0 > -1.0) else 1.0
                start = tr.t[k - 1] + w * (tr.t[k] - tr.t[k - 1])
            elif g0 >= 0 > g1 and start is not None:   # exit
                w = g0 / (g0 - g1) if g1 > -1.0 else 0.0
                end = tr.t[k - 1] + w * (tr.t[k] - tr.t[k - 1])
                if end - start > best:
                    best, best_start, best_end = end - start, start, end
                start = None
        if start is not None and g[-1] >= 0 and tr.t[-1] - start > best:
            best, best_start, best_end = tr.t[-1] - start, start, tr.t[-1]
        return {"duration": best, "t_start": best_start, "t_end": best_end,
                "met": best >= self.min_duration}


def window_state(tr, t):
    """Linear interpolation of the main flight variables at time t."""
    k = max(1, min(len(tr.t) - 1, next((i for i in range(len(tr.t)) if tr.t[i] >= t), len(tr.t) - 1)))
    t0, t1 = tr.t[k - 1], tr.t[k]
    w = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
    out = {}
    for f in ("alt", "mach", "q", "v", "fpa", "rho", "lat", "lon"):
        arr = getattr(tr, f)
        out[f] = arr[k - 1] + w * (arr[k] - arr[k - 1])
    return out


def state_at_altitude(tr, alt_km, phase="descent"):
    """Linear interpolation of the trajectory at a given altitude crossing."""
    ta, _, ia = tr.apogee()
    rng = range(ia, len(tr.t) - 1) if phase == "descent" else range(0, ia)
    h = alt_km * 1000.0
    for k in rng:
        a0, a1 = tr.alt[k], tr.alt[k + 1]
        if (a0 - h) * (a1 - h) <= 0 and a0 != a1:
            w = (h - a0) / (a1 - a0)
            out = {}
            for f in ("t", "v", "mach", "q", "fpa", "lat", "lon", "rho"):
                arr = getattr(tr, f)
                out[f] = arr[k] + w * (arr[k + 1] - arr[k])
            return out
    return None


def coverage(tr, stations, min_el=3.0, link=True, idx=None):
    """Elevation, slant range, visibility and link margin per station
    (for all samples, or only for the sample indices idx)."""
    res = {}
    ii = range(len(tr.t)) if idx is None else idx
    for name in stations:
        st = earth.STATIONS[name]
        el, rg, vis, lm = [], [], [], []
        for x, y, z in ((tr.x[k], tr.y[k], tr.z[k]) for k in ii):
            e, a, r = earth.look_angles(st, x, y, z)
            el.append(e)
            rg.append(r)
            vis.append(e >= min_el)
            lm.append(earth.link_margin_db(r) if link else float("nan"))
        res[name] = {"el": el, "range": rg, "visible": vis, "margin": lm}
    return res


def visible_fraction(tr, cov, name, mask=None):
    tot, vis = 0.0, 0.0
    for k in range(1, len(tr.t)):
        if mask is not None and not mask[k]:
            continue
        dt = tr.t[k] - tr.t[k - 1]
        tot += dt
        if cov[name]["visible"][k]:
            vis += dt
    return vis / tot if tot > 0 else float("nan")


def downrange_km(tr, site):
    return [earth.great_circle_km(site[0], site[1], la, lo) for la, lo in zip(tr.lat, tr.lon)]
