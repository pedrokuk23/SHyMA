"""Parametric geometry of an axisymmetric finned vehicle and its mass properties.

A configuration is
  * a body of revolution given by its meridian profile, a list of (x, r)
    points from the nose tip (x = 0) to the base, joined by straight segments
    (cones, frustums, cylinders); a curved nose is represented by several
    points;
  * any number of fin sets, each of n identical trapezoidal panels equally
    spaced in roll, defined by the exposed root chord, tip chord, exposed span,
    the axial offset of the tip leading edge (sweep) and the root leading-edge
    station;
  * the reference diameter (normally the largest body diameter).

The same objects describe every flight configuration of a multistage vehicle:
`Configuration.truncated(x_base, keep_fins)` returns the stack left after a
stage drop. `MassItem` lists give the centre of gravity and pitch inertia as
functions of time. A body profile can also be produced from a triangulated
surface (STL) by `profile_from_points`, which keeps the outer radius per axial
station (fins are then entered as fin sets).
"""
import math


class FinSet:
    def __init__(self, name, x_le, c_root, c_tip, span, x_tip_offset, n=4, t_c=0.03):
        self.name = name
        self.x_le = float(x_le)            # root leading edge, from the nose tip (m)
        self.c_root = float(c_root)        # exposed root chord (m)
        self.c_tip = float(c_tip)          # tip chord (m)
        self.span = float(span)            # exposed span of one panel (m)
        self.x_tip_offset = float(x_tip_offset)   # tip LE aft of root LE (m)
        self.n = int(n)
        self.t_c = float(t_c)

    @property
    def area(self):                      # one panel
        return 0.5 * (self.c_root + self.c_tip) * self.span

    @property
    def aspect(self):                    # two panels joined at the root (as in linear wing theory)
        return 2.0 * self.span ** 2 / self.area

    @property
    def taper(self):
        return self.c_tip / self.c_root

    def sweep_le(self):
        return math.atan2(self.x_tip_offset, self.span)

    def sweep_mid(self):
        return math.atan2(self.x_tip_offset + 0.5 * self.c_tip - 0.5 * self.c_root, self.span)

    def sweep_te(self):
        return math.atan2(self.x_tip_offset + self.c_tip - self.c_root, self.span)

    def mac(self):
        """Mean aerodynamic chord and the axial station of its leading edge."""
        lam = self.taper
        c = 2.0 / 3.0 * self.c_root * (1 + lam + lam * lam) / (1 + lam)
        y = self.span / 3.0 * (1 + 2 * lam) / (1 + lam)
        return c, self.x_le + y * self.x_tip_offset / self.span

    def centroid_x(self):
        """Axial station of the planform area centroid."""
        cr, ct, xt = self.c_root, self.c_tip, self.x_tip_offset
        # trapezoid with vertices (0,0) (cr,0) (xt+ct, s) (xt, s), x along the axis
        num = (cr * cr + cr * ct + ct * ct + xt * (cr + 2 * ct)) / 3.0
        return self.x_le + num / (cr + ct)

    def x_te(self):
        return self.x_le + self.c_root


class Configuration:
    def __init__(self, profile, fins, d_ref, name="vehicle"):
        self.profile = [(float(x), float(r)) for x, r in profile]
        self.fins = list(fins)
        self.d_ref = float(d_ref)
        self.name = name
        for (x0, _), (x1, _) in zip(self.profile, self.profile[1:]):
            if x1 < x0:
                raise ValueError("profile stations must increase")

    @property
    def s_ref(self):
        return math.pi * self.d_ref ** 2 / 4.0

    @property
    def length(self):
        return self.profile[-1][0] - self.profile[0][0]

    def radius(self, x):
        p = self.profile
        if x <= p[0][0]:
            return p[0][1]
        for (x0, r0), (x1, r1) in zip(p, p[1:]):
            if x <= x1:
                return r0 if x1 == x0 else r0 + (r1 - r0) * (x - x0) / (x1 - x0)
        return p[-1][1]

    def base_radius(self):
        return self.profile[-1][1]

    def segments(self):
        return [(x0, r0, x1, r1) for (x0, r0), (x1, r1) in zip(self.profile, self.profile[1:]) if x1 > x0]

    def volume(self):
        return sum(math.pi * (x1 - x0) * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0 for x0, r0, x1, r1 in self.segments())

    def planform_area(self):
        return sum((x1 - x0) * (r0 + r1) for x0, r0, x1, r1 in self.segments())

    def wetted_area(self):
        return sum(math.pi * (r0 + r1) * math.hypot(x1 - x0, r1 - r0) for x0, r0, x1, r1 in self.segments())

    def truncated(self, x_base, keep_fins, name=None):
        """Stack forward of station x_base (e.g. after the lower stage is dropped)."""
        pts = [(x, r) for x, r in self.profile if x < x_base]
        pts.append((x_base, self.radius(x_base)))
        fins = [f for f in self.fins if f.name in keep_fins]
        return Configuration(pts, fins, self.d_ref, name or self.name)


def profile_from_points(points, n_stations=200, axis=0):
    """Meridian profile (x, r) of a body of revolution from surface points
    (e.g. STL vertices): the axis is coordinate `axis`; r is the largest
    distance from the axis among the points of each axial bin. Fins appear as
    local radius excursions and should be removed (or entered as fin sets)."""
    xs = [p[axis] for p in points]
    x0, x1 = min(xs), max(xs)
    h = (x1 - x0) / n_stations
    rmax = [0.0] * (n_stations + 1)
    for p in points:
        k = min(int((p[axis] - x0) / h + 0.5), n_stations)
        r = math.sqrt(sum(p[i] ** 2 for i in range(3) if i != axis))
        rmax[k] = max(rmax[k], r)
    return [(x0 + k * h - x0, rmax[k]) for k in range(n_stations + 1)]


class MassItem:
    """Lumped mass of a vehicle component: a uniform rod of length `length`
    centred at station x (from the nose tip), present while `active(t)`."""

    def __init__(self, name, mass, x, length, active=None):
        self.name = name
        self.mass = mass                     # float or function of t
        self.x = x
        self.length = length
        self.active = active or (lambda t: True)

    def m(self, t, tc=None):
        if not self.active(t if tc is None else tc):
            return 0.0
        return self.mass(t) if callable(self.mass) else self.mass


def mass_properties(items, t, tc=None):
    """Total mass, centre of gravity (from the nose tip) and pitch inertia about the CG
    (stage drops decided at the control time tc, as in vehicle.Vehicle.mass)."""
    ms = [(it.m(t, tc), it) for it in items]
    m = sum(mi for mi, _ in ms)
    xc = sum(mi * it.x for mi, it in ms) / m
    I = sum(mi * (it.length ** 2 / 12.0 + (it.x - xc) ** 2) for mi, it in ms)
    return m, xc, I


# ---------------------------------------------------------------------------------------------
# VSB-30 estimate from public sources
# ---------------------------------------------------------------------------------------------
# Motor diameter 0.557 m and motor lengths: Palmerio et al. (2003), as tabulated by da Luz (2020).
# Axial layout, nose shape, payload diameter and fin planforms: 1:96 scale model of VSB-30 V05
# (N. J. Knudsen, 2015), measured from its vector outlines; axial stations scaled to the 12.6 m
# overall length; payload diameter 0.438 m (TEXUS standard; the model gives 0.43 m).
# Fin thickness and the number of fins per stage (3 in the model) are assumptions.
VSB30_LENGTH = 12.6
_MODEL_SECTIONS = [  # (name, model length m, diameter at the aft end m) from nose tip
    ("nose cap", 0.132, 0.133), ("nose", 0.725, 0.367), ("nose base", 0.325, 0.431),
    ("payload", 3.651, 0.431), ("adapter", 0.744, 0.517),
    ("S30 motor", 3.062, 0.517), ("S30 aft", 0.596, 0.517),
    ("S31 motor", 3.062, 0.517), ("S31 aft", 0.596, 0.517)]


def vsb30_geometry(n_fins=3, t_c=0.03, fin_scale=1.0):
    k = VSB30_LENGTH / sum(s[1] for s in _MODEL_SECTIONS)
    rad = {0.133: 0.133 * 0.438 / 0.431, 0.367: 0.367 * 0.438 / 0.431,
           0.431: 0.438, 0.517: 0.557}
    x, prof, st = 0.0, [(0.0, 0.090 * 0.438 / 0.431 / 2.0)], {}
    for name, L, d in _MODEL_SECTIONS:
        x += L * k
        prof.append((x, rad[d] / 2.0))
        st[name] = x
    s = 0.352778e-3 * 96.0 * fin_scale     # model points to full-scale metres
    s30 = FinSet("S30", st["S30 aft"] - 20.8 * s, 20.8 * s, 20.8 * s, 17.9 * s, 10.7 * s, n_fins, t_c)
    s31 = FinSet("S31", st["S31 aft"] - 19.2 * s, 19.2 * s, 11.9 * s, 13.1 * s, 7.3 * s, n_fins, t_c)
    full = Configuration(prof, [s30, s31], 0.557, "VSB-30")
    upper = full.truncated(st["S30 aft"], {"S30"}, "S30 + payload")
    return {"full": full, "upper": upper, "stations": st}


def vsb30_mass_items(veh, st):
    """Mass items of a SHyMA VSB-30 vehicle (vehicle.vsb30) on the geometry stations `st`:
    payload (nose and payload sections), and, for each motor, propellant (uniform over the motor
    case, burning at the stage's mass-flow rate) and dry mass (motor case, adapter, aft section
    and fins)."""
    s31, s30 = veh.stages
    x_pl0, x_pl1 = 0.0, st["payload"]
    items = [MassItem("payload", veh.payload_mass, 0.5 * (x_pl0 + x_pl1), x_pl1 - x_pl0)]

    def prop(stg):
        return lambda t: max(stg.prop_mass - stg.burned(min(max(t - stg.t_ign, 0.0), stg.burn_time)), 0.0)
    a30, b30 = st["adapter"], st["S30 motor"]
    items.append(MassItem("S30 propellant", prop(s30), 0.5 * (a30 + b30), b30 - a30,
                          active=lambda t, s=s30: t < s.t_drop))
    items.append(MassItem("S30 dry", s30.dry_mass, 0.5 * (st["payload"] + st["S30 aft"]) + 0.6,
                          st["S30 aft"] - st["payload"], active=lambda t, s=s30: t < s.t_drop))
    a31, b31 = st["S30 aft"], st["S31 motor"]
    items.append(MassItem("S31 propellant", prop(s31), 0.5 * (a31 + b31), b31 - a31,
                          active=lambda t, s=s31: t < s.t_drop))
    items.append(MassItem("S31 dry", s31.dry_mass, 0.5 * (st["S30 aft"] + st["S31 aft"]) + 0.6,
                          st["S31 aft"] - st["S30 aft"], active=lambda t, s=s31: t < s.t_drop))
    return items
