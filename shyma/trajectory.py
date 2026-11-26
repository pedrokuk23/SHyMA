"""3-DOF point-mass trajectory over a rotating Earth with J2 gravity.

State is integrated in an Earth-centred inertial frame (aligned with ECEF
at t = 0). Thrust acts along the air-relative velocity (unguided,
spin-stabilised, weathercock-stable rocket); drag opposes it. The
air-relative velocity accounts for Earth rotation and, optionally, for a
horizontal wind profile. A fixed-step RK4 scheme with a fine step inside the
sensible atmosphere and a coarse step in vacuum is used, with steps clamped
to land exactly on staging events and discrete events evaluated at the step
midpoint. The Earth model (spherical/WGS-84, J2, rotation) is selected with
earth.set_model().
"""
import math

from . import earth
from .atmosphere import atmosphere as _default_atmosphere


class Env:
    """Environment options for one propagation."""

    def __init__(self, level="mean", rho_factor=1.0, wind=None, atm=None, dT=None):
        self.dT = dT          # callable(alt_m) -> temperature offset (K), or None
        self.level = level
        self.rho_factor = rho_factor
        self.wind = wind      # callable(alt_m) -> (east, north) wind in m/s, or None
        self.atm = atm        # callable(alt_km, level, rho_factor) -> (T, p, rho, a), or None

    def atmosphere(self, alt_km):
        f = self.atm if self.atm is not None else _default_atmosphere
        T, p, rho, a = f(alt_km, self.level, self.rho_factor)
        if self.dT is not None:
            T2 = max(T + self.dT(alt_km * 1000.0), 100.0)
            pr = _hydrostatic_ratio(self.dT, f, self.level, self.rho_factor, alt_km)
            p *= pr                # pressure re-integrated with the perturbed column
            rho *= pr * T / T2
            a *= (T2 / T) ** 0.5
            T = T2
        return T, p, rho, a


_HYD_DZ = 0.25        # km
_HYD_TOP = 150.0      # km


def _hydrostatic_ratio(dT, f, level, rho_factor, alt_km):
    """Pressure ratio p'/p at alt_km when the temperature offset dT(h) is
    applied to the whole column above a fixed surface pressure:
    ln(p'/p) = int_0^h g dT / (R T (T + dT)) dh. The table is built once per
    offset object and base atmosphere and interpolated linearly."""
    key = (id(f), level, rho_factor)
    cache = getattr(dT, "_hyd_cache", None)
    if cache is None:
        cache = {}
        try:
            dT._hyd_cache = cache
        except AttributeError:
            pass
    tab = cache.get(key)
    if tab is None:
        R, g0, Re = 287.053, 9.80665, 6356.766
        n = int(_HYD_TOP / _HYD_DZ) + 1
        tab = [0.0] * n
        acc = 0.0
        prev = None
        for i in range(n):
            h = i * _HYD_DZ
            T = f(h, level, rho_factor)[0]
            T2 = max(T + dT(h * 1000.0), 100.0)
            g = g0 * (Re / (Re + h)) ** 2
            val = g * (T2 - T) / (R * T * T2)
            if prev is not None:
                acc += 0.5 * (prev + val) * _HYD_DZ * 1000.0
            prev = val
            tab[i] = acc
        cache[key] = tab
    x = min(max(alt_km, 0.0), _HYD_TOP) / _HYD_DZ
    i = min(int(x), len(tab) - 2)
    u = x - i
    return math.exp(tab[i] + u * (tab[i + 1] - tab[i]))


class Trajectory:
    """Container of sampled time histories (plain lists)."""

    FIELDS = ("t", "x", "y", "z", "alt", "lat", "lon", "v", "mach", "q", "rho",
              "fpa", "mass", "thrust", "acc", "vn", "ve", "vd")

    def __init__(self):
        for f in self.FIELDS:
            setattr(self, f, [])

    def append(self, **kw):
        for f in self.FIELDS:
            getattr(self, f).append(kw[f])

    def apogee(self):
        """(t, alt, index). The time and altitude are refined by a parabola
        through the three samples around the maximum (or taken from the
        exact per-step maximum tracked during integration, if available)."""
        i = max(range(len(self.alt)), key=lambda k: self.alt[k])
        ex = getattr(self, "apo_exact", None)
        if ex is not None and ex[1] >= self.alt[i]:
            return ex[0], ex[1], i
        if 0 < i < len(self.alt) - 1:
            t0, t1, t2 = self.t[i - 1], self.t[i], self.t[i + 1]
            a0, a1, a2 = self.alt[i - 1], self.alt[i], self.alt[i + 1]
            d = (t0 - t1) * (t0 - t2) * (t1 - t2)
            if d != 0:
                A = (t2 * (a1 - a0) + t1 * (a0 - a2) + t0 * (a2 - a1)) / d
                B = (t2 ** 2 * (a0 - a1) + t1 ** 2 * (a2 - a0) + t0 ** 2 * (a1 - a2)) / d
                C = (t1 * t2 * (t1 - t2) * a0 + t2 * t0 * (t2 - t0) * a1 + t0 * t1 * (t0 - t1) * a2) / d
                if A < 0:
                    tv = -B / (2 * A)
                    return tv, C - B * B / (4 * A), i
        return self.t[i], self.alt[i], i

    def impact(self):
        return self.lat[-1], self.lon[-1], self.t[-1]


def _air_velocity(t, x, y, z, vx, vy, vz, alt, env):
    """Air-relative velocity in ECI: v - omega x r - wind."""
    ux = vx + earth.OMEGA * y
    uy = vy - earth.OMEGA * x
    uz = vz
    if env.wind is not None and alt < 100e3:
        we, wn = env.wind(alt)
        if we or wn:
            xe, ye, ze = earth.eci_to_ecef(x, y, z, t)
            lat, lon, _ = earth.ecef_to_llh(xe, ye, ze)
            e, n, _u = earth.enu_basis(lat, lon)
            wxe = we * e[0] + wn * n[0]
            wye = we * e[1] + wn * n[1]
            wze = we * e[2] + wn * n[2]
            th = earth.OMEGA * t
            c, s = math.cos(th), math.sin(th)
            ux -= c * wxe - s * wye
            uy -= s * wxe + c * wye
            uz -= wze
    return ux, uy, uz


def _deriv(t, s, veh, env, tc=None):
    x, y, z, vx, vy, vz = s
    alt = earth.altitude(x, y, z)
    gx, gy, gz = earth.gravity(x, y, z)
    ux, uy, uz = _air_velocity(t, x, y, z, vx, vy, vz, alt, env)
    vr = math.sqrt(ux * ux + uy * uy + uz * uz)
    m = veh.mass(t, tc)
    ax, ay, az = gx, gy, gz
    p_amb = 0.0
    if alt < 300e3 and vr > 1e-6:
        T, p, rho, a = env.atmosphere(alt / 1000.0)
        p_amb = p
        if rho > 0.0:
            cda = (veh.drag_area(t, vr / a, tc, rho * vr * (T + 110.4) / (1.458e-6 * T ** 1.5))
                   if getattr(veh, "cd_model", None) is not None else veh.drag_area(t, vr / a, tc))
            k = -0.5 * rho * cda * vr / m
            ax += k * ux
            ay += k * uy
            az += k * uz
    F = veh.thrust(t, tc, p_amb)
    if F > 0.0 and vr > 1e-6:
        k = F / (m * vr)
        ax += k * ux
        ay += k * uy
        az += k * uz
    return (vx, vy, vz, ax, ay, az)


def _rk4(t, s, h, veh, env):
    tc = t + 0.5 * h
    k1 = _deriv(t, s, veh, env, tc)
    s2 = [s[i] + 0.5 * h * k1[i] for i in range(6)]
    k2 = _deriv(t + 0.5 * h, s2, veh, env, tc)
    s3 = [s[i] + 0.5 * h * k2[i] for i in range(6)]
    k3 = _deriv(t + 0.5 * h, s3, veh, env, tc)
    s4 = [s[i] + h * k3[i] for i in range(6)]
    k4 = _deriv(t + h, s4, veh, env, tc)
    return [s[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6)]


def _record(tr, t, s, veh, env):
    x, y, z, vx, vy, vz = s
    alt = earth.altitude(x, y, z)
    ux, uy, uz = _air_velocity(t, x, y, z, vx, vy, vz, alt, env)
    vr = math.sqrt(ux * ux + uy * uy + uz * uz)
    T, p, rho, a = env.atmosphere(max(alt, 0.0) / 1000.0)
    xe, ye, ze = earth.eci_to_ecef(x, y, z, t)
    lat, lon, _ = earth.ecef_to_llh(xe, ye, ze)
    e, n, u = earth.enu_basis(lat, lon)
    # Earth-relative velocity (no wind) in local NED
    gx, gy = vx + earth.OMEGA * y, vy - earth.OMEGA * x
    th = earth.OMEGA * t
    c, sn = math.cos(th), math.sin(th)
    vxe, vye, vze = c * gx + sn * gy, -sn * gx + c * gy, vz
    ve_ = vxe * e[0] + vye * e[1] + vze * e[2]
    vn_ = vxe * n[0] + vye * n[1] + vze * n[2]
    vu_ = vxe * u[0] + vye * u[1] + vze * u[2]
    vg = math.sqrt(ve_ * ve_ + vn_ * vn_ + vu_ * vu_)
    fpa = math.degrees(math.asin(max(-1.0, min(1.0, vu_ / vg)))) if vg > 0 else 90.0
    d = _deriv(t, s, veh, env)
    g = earth.gravity(x, y, z)
    acc = math.sqrt(sum((d[3 + i] - g[i]) ** 2 for i in range(3)))
    tr.append(t=t, x=xe, y=ye, z=ze, alt=alt, lat=lat, lon=lon, v=vr,
              mach=vr / a if a > 0 else 0.0, q=0.5 * rho * vr * vr, rho=rho, fpa=fpa,
              mass=veh.mass(t), thrust=veh.thrust(t), acc=acc, vn=vn_, ve=ve_, vd=-vu_)


def launch_state(site, elevation_deg, azimuth_deg, veh, rail_length=12.0):
    """State at rail exit for a rail-launched rocket (analytic rail phase).

    Returns (t_exit, state) with the rail exit velocity along the rail and the
    Earth-rotation velocity added. Mass and thrust variation on the rail are
    neglected (rail time < 1 s); the thrust includes the ambient-pressure loss
    at the launcher altitude, as in the integration."""
    lat, lon, h = site
    m0 = veh.mass(0.0)
    F = veh.thrust(0.0, p_amb=_default_atmosphere(max(h, 0.0) / 1000.0)[1])
    acc = F / m0 - 9.80665 * math.sin(math.radians(elevation_deg))
    v_exit = math.sqrt(2.0 * acc * rail_length)
    t_exit = v_exit / acc
    e, n, u = earth.enu_basis(lat, lon)
    el, az = math.radians(elevation_deg), math.radians(azimuth_deg)
    d = [math.cos(el) * math.sin(az) * e[i] + math.cos(el) * math.cos(az) * n[i] + math.sin(el) * u[i]
         for i in range(3)]
    x, y, z = earth.llh_to_ecef(lat, lon, h)
    x, y, z = x + d[0] * rail_length, y + d[1] * rail_length, z + d[2] * rail_length
    return t_exit, ecef_to_eci_state(x, y, z, [v_exit * di for di in d], t_exit)


def ecef_to_eci_state(x, y, z, v_rel_ecef, t):
    """Inertial state from an ECEF position and Earth-relative ECEF velocity."""
    th = earth.OMEGA * t
    c, s = math.cos(th), math.sin(th)
    xi, yi = c * x - s * y, s * x + c * y
    vxi = c * v_rel_ecef[0] - s * v_rel_ecef[1] - earth.OMEGA * yi
    vyi = s * v_rel_ecef[0] + c * v_rel_ecef[1] + earth.OMEGA * xi
    return [xi, yi, z, vxi, vyi, v_rel_ecef[2]]


def state_from_ned(lat, lon, h, vn, ve, vd, t=0.0):
    """Inertial state from geodetic position and Earth-relative NED velocity."""
    x, y, z = earth.llh_to_ecef(lat, lon, h)
    e, n, u = earth.enu_basis(lat, lon)
    v = [vn * n[i] + ve * e[i] - vd * u[i] for i in range(3)]
    return ecef_to_eci_state(x, y, z, v, t)


def state_from_apogee(site_lat, site_lon, alt_m, v_horiz, heading_deg, t0=0.0):
    """Inertial state at an apogee point with a given horizontal (relative)
    speed and heading. Used to reconstruct a re-entry from flight data."""
    hd = math.radians(heading_deg)
    return state_from_ned(site_lat, site_lon, alt_m, v_horiz * math.cos(hd),
                          v_horiz * math.sin(hd), 0.0, t0)


def propagate(veh, t0, s0, level="mean", rho_factor=1.0, dt_atm=0.02, dt_vac=0.5,
              t_max=3000.0, record_every=1, stop_alt=0.0, stop_after_apogee_alt=None,
              wind=None, atm=None, t_end=None, dT=None):
    """Integrate until ground impact (alt <= stop_alt after lift-off) or t_end."""
    env = Env(level, rho_factor, wind, atm, dT)
    tr = Trajectory()
    events = [e for e in veh.event_times() if e > t0]
    if t_end is not None:
        events = sorted(set(events + [t_end]))
        t_max = t_end
    t, s = t0, list(s0)
    _record(tr, t, s, veh, env)
    n = 0
    climbing = True
    prev_alt = earth.altitude(*s[:3])
    apo = (t, prev_alt)
    s_pp, t_pp = None, None
    while t < t_max - 1e-9:
        alt = earth.altitude(*s[:3])
        h = dt_atm if (alt < 110e3 or veh.thrust(t) > 0.0) else dt_vac
        while events and events[0] <= t + 1e-9:
            events.pop(0)
        if events and t + h > events[0]:
            h = events[0] - t
        s_prev, t_prev = s, t
        s = _rk4(t, s, h, veh, env)
        t += h
        n += 1
        alt = earth.altitude(*s[:3])
        if climbing and alt < prev_alt:
            climbing = False
            # apogee between t_prev and t: cubic Hermite on the radial motion
            apo = _apogee_hermite(t_prev, s_prev, t, s, apo)
            if s_pp is not None:   # the maximum may lie in the previous step
                apo = _apogee_hermite(t_pp, s_pp, t_prev, s_prev, apo)
        if alt <= stop_alt and t > t0 + 1.0:
            # impact (or stop altitude): linear interpolation within the step
            w = (prev_alt - stop_alt) / (prev_alt - alt) if prev_alt != alt else 1.0
            s = [a + w * (b - a) for a, b in zip(s_prev, s)]
            t = t_prev + w * (t - t_prev)
            _record(tr, t, s, veh, env)
            break
        prev_alt = alt
        s_pp, t_pp = s_prev, t_prev
        if n % record_every == 0 or t >= t_max - 1e-9:
            _record(tr, t, s, veh, env)
        if stop_after_apogee_alt is not None and not climbing and alt <= stop_after_apogee_alt:
            _record(tr, t, s, veh, env)
            break
    tr.final = (t, s)
    if not climbing:
        tr.apo_exact = apo
    return tr


def _apogee_hermite(t0, s0, t1, s1, default):
    """Time and altitude of the altitude maximum within one step, from a
    cubic Hermite interpolant of position (velocities as end slopes)."""
    best = (t0, earth.altitude(*s0[:3]))
    h = t1 - t0
    for k in range(1, 41):
        u = k / 40.0
        h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
        h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
        p = [h00 * s0[i] + h10 * h * s0[i + 3] + h01 * s1[i] + h11 * h * s1[i + 3] for i in range(3)]
        a = earth.altitude(*p)
        if a > best[1]:
            best = (t0 + u * h, a)
    return best if best[1] >= default[1] else default
