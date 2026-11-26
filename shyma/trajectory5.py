"""Pseudo-5-DOF propagation: point-mass translation plus pitch/yaw rotation of
the body axis (roll/spin neglected), for unguided, statically stable rockets.

The 3-DOF model of trajectory.py assumes instantaneous weathercocking (thrust
always along the air-relative velocity). At low speed this overstates the
turn of the rocket into the wind and the gravity turn after rail exit. Here
the body axis b is a state: it starts along the launcher and rotates under
the aerodynamic restoring moment of the normal force acting at the centre of
pressure, with pitch damping. Thrust acts along b.

State: [r(3), v(3), b(3), w(3)] in ECI; w is the transverse angular velocity.
"""
import math

from . import earth, trajectory
from .trajectory import Env, Trajectory, _air_velocity


class Aero5:
    """Generic stability data of a finned slender rocket.

    cn_alpha  : normal-force slope (1/rad, referred to ref_area)
    margin_cal: static margin in calibers (CP aft of CG)
    cmq       : pitch-damping coefficient (1/rad, referred to ref_area and diameter)
    """

    def __init__(self, diameter, cn_alpha=12.0, margin_cal=1.5, cmq=-300.0,
                 lengths=(12.6, 9.4, 4.0)):
        self.d = diameter
        self.cn_alpha = cn_alpha
        self.margin = margin_cal * diameter
        self.cmq = cmq
        self.lengths = lengths   # stack length: full, after first-stage drop, payload alone

    def length(self, veh, t, tc=None):
        tc = t if tc is None else tc
        dropped = sum(1 for s in veh.stages if tc >= s.t_drop)
        return self.lengths[min(dropped, len(self.lengths) - 1)]

    def inertia(self, veh, t, tc=None):
        L = self.length(veh, t, tc)
        return veh.mass(t, tc) * L * L / 12.0


def _norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return (v[0] / n, v[1] / n, v[2] / n), n


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _deriv5(t, s, veh, aero, env, tc=None):
    x, y, z, vx, vy, vz = s[0:6]
    b = s[6:9]
    w = s[9:12]
    alt = earth.altitude(x, y, z)
    gx, gy, gz = earth.gravity(x, y, z)
    ux, uy, uz = _air_velocity(t, x, y, z, vx, vy, vz, alt, env)
    vr = math.sqrt(ux * ux + uy * uy + uz * uz)
    m = veh.mass(t, tc)
    ax, ay, az = gx, gy, gz
    tau = (0.0, 0.0, 0.0)
    p_amb = 0.0
    if alt < 300e3 and vr > 1e-3:
        T, p, rho, a = env.atmosphere(alt / 1000.0)
        p_amb = p
        if rho > 0.0:
            q = 0.5 * rho * vr * vr
            uh = (ux / vr, uy / vr, uz / vr)
            # axial drag along -uh
            cda = (veh.drag_area(t, vr / a, tc, rho * vr * (T + 110.4) / (1.458e-6 * T ** 1.5))
                   if getattr(veh, "cd_model", None) is not None else veh.drag_area(t, vr / a, tc))
            k = -q * cda / m
            ax += k * uh[0]
            ay += k * uh[1]
            az += k * uh[2]
            # crossflow component of the flight direction, perpendicular to b
            ub = uh[0] * b[0] + uh[1] * b[1] + uh[2] * b[2]
            up = (uh[0] - ub * b[0], uh[1] - ub * b[1], uh[2] - ub * b[2])
            sa = math.sqrt(up[0] ** 2 + up[1] ** 2 + up[2] ** 2)   # sin(alpha)
            S = veh.ref_area
            if hasattr(aero, "coeffs"):      # geometry-based, Mach- and time-dependent derivatives
                cna, margin, cmq = aero.coeffs(veh, t, vr / a, tc)
            else:
                cna, margin, cmq = aero.cn_alpha, aero.margin, aero.cmq
            if sa > 1e-12:
                fn = q * S * cna * sa
                upn = (up[0] / sa, up[1] / sa, up[2] / sa)
                ax -= fn / m * upn[0]
                ay -= fn / m * upn[1]
                az -= fn / m * upn[2]
                # restoring moment: rotates b towards the flight direction
                c = _cross(b, upn)
                mt = margin * fn
                tau = (mt * c[0], mt * c[1], mt * c[2])
            # pitch damping
            kd = q * S * aero.d * aero.d * cmq / (2.0 * vr)
            tau = (tau[0] + kd * w[0], tau[1] + kd * w[1], tau[2] + kd * w[2])
    F = veh.thrust(t, tc, p_amb)
    if F > 0.0:
        ax += F / m * b[0]
        ay += F / m * b[1]
        az += F / m * b[2]
    I = aero.inertia(veh, t, tc)
    dw = (tau[0] / I, tau[1] / I, tau[2] / I)
    db = _cross(w, b)
    return (vx, vy, vz, ax, ay, az, db[0], db[1], db[2], dw[0], dw[1], dw[2])


def _tilt(s, angle_deg, phase):
    """Tilt the body axis s[6:9] by angle_deg; phase 0 = in the vertical plane, upwards."""
    b = s[6:9]
    r, _ = _norm(s[0:3])
    rb = sum(r[i] * b[i] for i in range(3))
    e1, n1 = _norm((r[0] - rb * b[0], r[1] - rb * b[1], r[2] - rb * b[2]))
    e2 = _cross(b, e1)
    e = [math.cos(phase) * e1[i] + math.sin(phase) * e2[i] for i in range(3)]
    th = math.radians(angle_deg)
    nb, _ = _norm(tuple(math.cos(th) * b[i] + math.sin(th) * e[i] for i in range(3)))
    out = list(s)
    out[6:9] = nb
    return out


def _rk4(t, s, h, veh, aero, env):
    tc = t + 0.5 * h
    k1 = _deriv5(t, s, veh, aero, env, tc)
    s2 = [s[i] + 0.5 * h * k1[i] for i in range(12)]
    k2 = _deriv5(t + 0.5 * h, s2, veh, aero, env, tc)
    s3 = [s[i] + 0.5 * h * k2[i] for i in range(12)]
    k3 = _deriv5(t + 0.5 * h, s3, veh, aero, env, tc)
    s4 = [s[i] + h * k3[i] for i in range(12)]
    k4 = _deriv5(t + h, s4, veh, aero, env, tc)
    out = [s[i] + h / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(12)]
    bn, _ = _norm(out[6:9])
    out[6:9] = bn
    # keep w transverse to b (roll neglected)
    wb = sum(out[9 + i] * bn[i] for i in range(3))
    out[9:12] = [out[9 + i] - wb * bn[i] for i in range(3)]
    return out


def launch_state5(site, elevation_deg, azimuth_deg, veh, rail_length=12.0):
    t_exit, s = trajectory.launch_state(site, elevation_deg, azimuth_deg, veh, rail_length)
    # body axis along the rail (Earth-fixed direction), expressed in ECI at t_exit
    e, n, u = earth.enu_basis(site[0], site[1])
    el, az = math.radians(elevation_deg), math.radians(azimuth_deg)
    d = [math.cos(el) * math.sin(az) * e[i] + math.cos(el) * math.cos(az) * n[i] + math.sin(el) * u[i]
         for i in range(3)]
    th = earth.OMEGA * t_exit
    c, sn = math.cos(th), math.sin(th)
    b = [c * d[0] - sn * d[1], sn * d[0] + c * d[1], d[2]]
    # rail-imposed angular velocity: Earth rotation (transverse part)
    om = (0.0, 0.0, earth.OMEGA)
    ob = sum(om[i] * b[i] for i in range(3))
    w = [om[i] - ob * b[i] for i in range(3)]
    return t_exit, s + b + w


def propagate5(veh, aero, t0, s0, level="mean", rho_factor=1.0, dt_atm=0.01, dt_vac=0.5,
               t_max=3000.0, record_every=1, stop_alt=0.0, wind=None, atm=None,
               stop_after_apogee_alt=None, t_end=None, dT=None, att_kick=None):
    """att_kick = (t_kick, angle_deg, phase_rad): optional attitude error applied to the
    body axis at t_kick (e.g., upper-stage ignition after a long coast): the axis is
    tilted by angle_deg towards a direction at phase_rad from the local vertical plane."""
    env = Env(level, rho_factor, wind, atm, dT)
    tr = Trajectory()
    tr.alpha = []
    events = [e for e in veh.event_times() if e > t0]
    kick_done = att_kick is None
    if not kick_done:
        events = sorted(set(events + [att_kick[0]]))
    if t_end is not None:
        events = sorted(set(events + [t_end]))
        t_max = t_end
    t, s = t0, list(s0)

    def rec():
        trajectory._record(tr, t, s[:6], veh, env)
        d5 = _deriv5(t, s, veh, aero, env)          # replace the 3-DOF acceleration by the 5-DOF one
        g5 = earth.gravity(*s[:3])
        tr.acc[-1] = math.sqrt(sum((d5[3 + i] - g5[i]) ** 2 for i in range(3)))
        x, y, z, vx, vy, vz = s[:6]
        alt = earth.altitude(x, y, z)
        ux, uy, uz = _air_velocity(t, x, y, z, vx, vy, vz, alt, env)
        uh, vr = _norm((ux, uy, uz))
        cb = sum(uh[i] * s[6 + i] for i in range(3))
        tr.alpha.append(math.degrees(math.acos(max(-1.0, min(1.0, cb)))))
    rec()
    n = 0
    climbing = True
    prev_alt = earth.altitude(*s[:3])
    while t < t_max - 1e-9:
        alt = earth.altitude(*s[:3])
        h = dt_atm if (alt < 110e3 or veh.thrust(t) > 0.0) else dt_vac
        while events and events[0] <= t + 1e-9:
            events.pop(0)
        if events and t + h > events[0]:
            h = events[0] - t
        s = _rk4(t, s, h, veh, aero, env)
        t += h
        n += 1
        if not kick_done and t >= att_kick[0] - 1e-9:
            s = _tilt(s, att_kick[1], att_kick[2])
            kick_done = True
        alt = earth.altitude(*s[:3])
        if climbing and alt < prev_alt:
            climbing = False
        prev_alt = alt
        if n % record_every == 0 or alt <= stop_alt or t >= t_max - 1e-9:
            rec()
        if alt <= stop_alt and t > t0 + 1.0:
            break
        if stop_after_apogee_alt is not None and not climbing and alt <= stop_after_apogee_alt:
            rec()
            break
    tr.final = (t, s)
    return tr
