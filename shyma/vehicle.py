"""Vehicle definition for SHyMA: unguided multistage sounding rocket
with constant-average-thrust solid motors and a generic drag curve."""
import bisect
import math

G0 = 9.80665

# Generic slender sounding-rocket zero-lift drag curve (power-off), used as a
# conceptual-phase default. It is scaled by Vehicle.cd_scale and per phase.
_CD_M = [0.0, 0.6, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 25.0]
_CD_V = [0.32, 0.32, 0.38, 0.50, 0.55, 0.54, 0.48, 0.40, 0.32, 0.27, 0.24, 0.22, 0.20, 0.19, 0.19]


def cd_generic(mach):
    if mach <= _CD_M[0]:
        return _CD_V[0]
    if mach >= _CD_M[-1]:
        return _CD_V[-1]
    k = bisect.bisect_right(_CD_M, mach)
    m0, m1 = _CD_M[k - 1], _CD_M[k]
    return _CD_V[k - 1] + (_CD_V[k] - _CD_V[k - 1]) * (mach - m0) / (m1 - m0)


class Stage:
    def __init__(self, name, prop_mass, dry_mass, burn_time, isp, t_ign, t_drop, shape=0.0, exit_area=0.0):
        self.name = name
        self.exit_area = exit_area   # nozzle exit area (m^2) for the ambient-pressure thrust loss; 0 = none
        self.shape = shape      # linear thrust-shape factor k: F = F_mean (1 + k (1 - 2 tau / t_b))
        self.prop_mass = prop_mass
        self.dry_mass = dry_mass
        self.burn_time = burn_time
        self.isp = isp
        self.t_ign = t_ign
        self.t_drop = t_drop

    @property
    def mdot(self):
        return self.prop_mass / self.burn_time

    @property
    def thrust(self):
        return self.mdot * self.isp * G0

    def burning(self, t):
        return self.t_ign <= t < self.t_ign + self.burn_time

    def thrust_at(self, t):
        """Thrust at time t (burning assumed). k > 0: regressive, k < 0:
        progressive; the total impulse and the burn time are unchanged."""
        if self.shape == 0.0:
            return self.thrust
        tau = min(max(t - self.t_ign, 0.0), self.burn_time)
        return self.thrust * (1.0 + self.shape * (1.0 - 2.0 * tau / self.burn_time))

    def burned(self, tau):
        """Propellant mass burned after tau seconds of burn."""
        if self.shape == 0.0:
            return self.mdot * tau
        return self.mdot * (tau + self.shape * (tau - tau * tau / self.burn_time))


class Vehicle:
    """Stack of stages plus payload. All stage dry masses stay attached until
    their drop time; the payload is never dropped."""

    def __init__(self, stages, payload_mass, ref_area, cd_scale=1.0,
                 payload_cd_scale=1.0, payload_area=None, name="vehicle"):
        self.stages = stages
        self.payload_mass = payload_mass
        self.ref_area = ref_area
        self.cd_scale = cd_scale
        self.payload_cd_scale = payload_cd_scale
        self.payload_area = payload_area if payload_area else ref_area
        self.name = name

    def event_times(self):
        ev = set()
        for s in self.stages:
            ev.update([s.t_ign, s.t_ign + s.burn_time, s.t_drop])
        return sorted(ev)

    def mass(self, t, tc=None):
        """Mass at time t. Discrete events (stage drops) are decided at the
        control time tc (the RK4 step midpoint), so that no integration step
        straddles a discontinuity."""
        tc = t if tc is None else tc
        m = self.payload_mass
        for s in self.stages:
            if tc >= s.t_drop:
                continue
            m += s.dry_mass
            tb = min(max(t, s.t_ign), s.t_ign + s.burn_time)
            if tc >= s.t_ign + s.burn_time:
                continue
            m += s.prop_mass - s.burned(tb - s.t_ign)
        return m

    def thrust(self, t, tc=None, p_amb=0.0):
        """Thrust (N); with p_amb > 0 and a stage exit area, the ambient-pressure
        loss p_amb A_e is subtracted (the specific impulse is then a vacuum value)."""
        tc = t if tc is None else tc
        return sum(max(s.thrust_at(t) - p_amb * s.exit_area, 0.0)
                   for s in self.stages if s.burning(tc) and tc < s.t_drop)

    def attached(self, t):
        return any(t < s.t_drop for s in self.stages)

    cd_model = None   # optional geometric drag of the attached stack, f(t, mach, tc, re_per_m) (shyma.drag)

    tv_excess = 0.0   # excess forebody drag of a blunt test vehicle while attached, in units of cd_generic*ref_area

    def drag_area(self, t, mach, tc=None, re_per_m=None):
        """Drag area C_D S (m^2). While the payload is attached, the stack drag is the generic
        (or geometric) rocket drag plus the excess drag of a test vehicle blunter than the
        generic slender nose (tv_excess), so that the stack never has less drag than the
        separated test vehicle carried on its nose."""
        tc_ = t if tc is None else tc
        if self.attached(tc_):
            extra = self.tv_excess * cd_generic(mach) * self.ref_area
            if self.cd_model is not None:
                return self.cd_model(t, mach, tc_, re_per_m) * self.cd_scale * self.ref_area + extra
            return cd_generic(mach) * self.cd_scale * self.ref_area + extra
        return cd_generic(mach) * self.payload_cd_scale * self.payload_area

    def total_time_burn_end(self):
        return max(s.t_ign + s.burn_time for s in self.stages)


def vsb30(payload_mass=400.0, isp=None, s31_burn=None, coast=2.0, t_payload_sep=60.0,
          cd_scale=1.0, impulse_factor=(1.0, 1.0), dry_factor=1.0,
          payload_cd_scale=1.0, beta_tv=None, burn_factor=(1.0, 1.0), exit_area=None,
          beta_stack=None):
    """VSB-30 (S31 + S30) model from public data.

    Stage gross/empty masses: 900/284 kg and 1200/341 kg; S30 action time
    32 s (Garcia et al., JATM 2011); diameter 0.557 m. Two parameters are
    calibrated on flight V01 (calibration.py): the effective specific impulse
    (common to both motors) and the effective burn time of the S31 booster,
    which represents its boost-dominated thrust curve with an equivalent
    constant thrust (the published 16 s action time includes the tail-off).

    beta_tv     : if given, ballistic coefficient m/(C_D S) of the separated
                  test vehicle at Mach 6 (kg/m^2); the generic C_D(M) shape is
                  kept and rescaled. If None, the payload uses the generic
                  rocket curve (as for the VSB-30 validation payloads). A test
                  vehicle blunter than the generic nose (C_D S above the generic
                  value on the 0.557 m reference) adds its excess drag to the
                  attached stack (Vehicle.tv_excess), also in captive concepts.
    beta_stack  : design coefficient used for that stack excess (default: beta_tv).
    burn_factor : multiplicative dispersion of the S31/S30 burn times at
                  constant total impulse (propellant-temperature effect).
    t_payload_sep = inf keeps the spent S30 attached (captive configuration).
    exit_area   : nozzle exit areas (S31, S30) in m^2 for the ambient-pressure
                  thrust loss; None = module EXIT_AREA (exit diameter EXIT_DIAMETER_M for both motors).
    """
    ae = EXIT_AREA if exit_area is None else exit_area
    if isp is None or s31_burn is None:
        if not CALIBRATED:
            raise RuntimeError("calibration.json missing: run python3 -m shyma.calibration")
    if isp is None:
        isp = ISP_VSB30
    if s31_burn is None:
        s31_burn = S31_BURN
    b1 = s31_burn * burn_factor[0]
    s31 = Stage("S31", 616.0 * impulse_factor[0], 284.0 * dry_factor, b1, isp,
                t_ign=0.0, t_drop=b1, shape=THRUST_SHAPE[0], exit_area=ae[0])
    t2 = b1 + coast
    s30 = Stage("S30", 859.0 * impulse_factor[1], 341.0 * dry_factor, S30_BURN * burn_factor[1], isp,
                t_ign=t2, t_drop=t_payload_sep, shape=THRUST_SHAPE[1], exit_area=ae[1])
    area = math.pi * 0.557 ** 2 / 4.0
    tv_excess = 0.0
    if beta_tv is not None:
        tv_scale = payload_mass / (beta_tv * area * cd_generic(6.0))
        payload_cd_scale = payload_cd_scale * tv_scale
        # the excess drag carried by the attached stack follows the design
        # coefficient (beta_stack) when given, so that a dispersed test-vehicle
        # coefficient changes only the free flight; the stack drag is
        # dispersed separately through the axial-drag scale
        b_st = beta_stack if beta_stack is not None else beta_tv
        tv_excess = max(payload_mass / (b_st * area * cd_generic(6.0)) - 1.0, 0.0)
    v = Vehicle([s31, s30], payload_mass, area, cd_scale=cd_scale,
                payload_cd_scale=payload_cd_scale, name="VSB-30")
    v.tv_excess = tv_excess
    v.t_burnout = t2 + S30_BURN * burn_factor[1]
    if DRAG_FACTORY is not None:          # optional geometric drag (shyma.drag); None = generic curve
        DRAG_FACTORY(v)
    return v


# Defaults; calibration.py fits them on flight V01 (236.0 km apogee, 189.0 km range, 392 kg, 84 deg);
# overwritten at import time if calibration.json exists.
ISP_VSB30 = 255.0
# linear thrust-shape factors (S31, S30); (0, 0) is the constant equivalent thrust of the paper
THRUST_SHAPE = (0.0, 0.0)
# Nozzle exit areas (S31, S30) for the ambient-pressure thrust loss F = mdot Isp_vac g0 - p_a A_e.
# No public nozzle data: exit diameter 0.45 m for both motors (0.40-0.50 m in the sensitivity study,
# about 0.8 of the 0.557 m motor diameter); the calibrated Isp is then a vacuum value.
EXIT_DIAMETER_M = 0.45
EXIT_AREA = (math.pi * EXIT_DIAMETER_M ** 2 / 4.0,) * 2
# optional callable veh -> None that attaches a geometric drag model (e.g.
# lambda v: shyma.drag.attach(v, shyma.drag.vsb30_drag())); None keeps the generic curve
DRAG_FACTORY = None
S31_BURN = 16.0
S30_BURN = 32.0
CALIBRATED = False
try:
    import json
    import os
    _p = os.path.join(os.path.dirname(__file__), "calibration.json")
    if os.path.exists(_p):
        with open(_p) as _f:
            _c = json.load(_f)
            ISP_VSB30 = _c["isp_vsb30"]
            S31_BURN = _c.get("s31_burn", 16.0)
            CALIBRATED = True
            _ae = _c.get("exit_area_m2")
            if _ae is not None and any(abs(a - b) > 1e-9 for a, b in zip(_ae, EXIT_AREA)):
                import warnings
                warnings.warn("calibration.json was produced with exit areas %s, not the current EXIT_AREA %s: "
                              "run python3 -m shyma.calibration" % (_ae, EXIT_AREA))
except Exception:  # pragma: no cover
    pass


class BallisticBody:
    """Unpowered body defined by its ballistic coefficient beta = m/(Cd A)."""

    def __init__(self, mass, beta, name="body"):
        self.m = mass
        self.beta = beta
        self.name = name
        self.stages = []

    def event_times(self):
        return []

    def mass(self, t, tc=None):
        return self.m

    def thrust(self, t, tc=None, p_amb=0.0):
        return 0.0

    def attached(self, t):
        return False

    def drag_area(self, t, mach, tc=None):
        return self.m / self.beta

    def total_time_burn_end(self):
        return 0.0
