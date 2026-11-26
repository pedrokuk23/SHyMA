"""High-level flight driver for VSB-30-boosted experiments.

Ascent (rail exit to payload separation) is propagated with the pseudo-5-DOF
model (body-axis dynamics, weathercocking); after separation the payload is
propagated as a point mass (3-DOF), assuming a nose-first, zero-angle-of-attack
attitude (actively or aerodynamically stabilised test vehicle).
"""
from . import trajectory, trajectory5, vehicle

VSB30_AERO = dict(cn_alpha=12.0, margin_cal=1.5, cmq=-300.0, lengths=(12.6, 9.4, 4.0))
# Optional factory veh -> geometry-based derivatives (e.g. shyma.stability.vsb30_aero); None keeps
# the generic constants of VSB30_AERO (the paper's baseline).
AERO_FACTORY = None


def merge(a, b):
    out = trajectory.Trajectory()
    for f in trajectory.Trajectory.FIELDS:
        setattr(out, f, getattr(a, f) + getattr(b, f)[1:])
    out.alpha = getattr(a, "alpha", []) + [0.0] * (len(b.t) - 1)
    out.final = b.final
    if hasattr(b, "apo_exact"):
        out.apo_exact = b.apo_exact
    return out


def fly(veh, site, el, az, t_sep, wind=None, atm=None, rho_factor=1.0, dt=0.01,
        record_every=2, aero=None, rail_length=12.0, model="5dof", stop_after_apogee_alt=None,
        dT=None, att_kick=None):
    if model == "3dof":
        te, s0 = trajectory.launch_state(site, el, az, veh, rail_length)
        return trajectory.propagate(veh, te, s0, rho_factor=rho_factor, dt_atm=dt,
                                    record_every=record_every, wind=wind, atm=atm, dT=dT,
                                    stop_after_apogee_alt=stop_after_apogee_alt)
    if aero is None and AERO_FACTORY is not None:
        aero = AERO_FACTORY(veh)
    if hasattr(aero, "coeffs"):          # geometry-based derivatives (shyma.stability)
        a = aero
    else:
        a = trajectory5.Aero5(veh.diameter if hasattr(veh, "diameter") else 0.557,
                              **(aero or VSB30_AERO))
    te, s0 = trajectory5.launch_state5(site, el, az, veh, rail_length)
    if att_kick is not None and att_kick[0] is None:     # attitude error just after rail exit (tip-off)
        att_kick = (te + 0.5, att_kick[1], att_kick[2])
    tr1 = trajectory5.propagate5(veh, a, te, s0, rho_factor=rho_factor, dt_atm=dt,
                                 record_every=record_every, wind=wind, atm=atm, t_end=t_sep, dT=dT,
                                 att_kick=att_kick)
    t1, s1 = tr1.final
    tr2 = trajectory.propagate(veh, t1, s1[:6], rho_factor=rho_factor, dt_atm=dt,
                               record_every=record_every, wind=wind, atm=atm, dT=dT,
                               stop_after_apogee_alt=stop_after_apogee_alt)
    return merge(tr1, tr2)


def state_at(veh, site, el, az, t, aero=None, rail_length=12.0, dt=0.01, wind=None):
    """ECI state (t, [x, y, z, vx, vy, vz]) of the ascending stack at time t
    (5-DOF ascent), e.g. at a stage drop."""
    if aero is None and AERO_FACTORY is not None:
        aero = AERO_FACTORY(veh)
    a = aero if hasattr(aero, "coeffs") else trajectory5.Aero5(0.557, **(aero or VSB30_AERO))
    te, s0 = trajectory5.launch_state5(site, el, az, veh, rail_length)
    tr = trajectory5.propagate5(veh, a, te, s0, dt_atm=dt, record_every=1000, wind=wind, t_end=t)
    return tr.final[0], tr.final[1][:6]


def spent_stage_impact(veh, site, el, az, t_drop, beta, wind=None):
    """Impact point of a spent stage dropped at t_drop, modelled as a
    tumbling ballistic body of ballistic coefficient beta (kg/m^2)."""
    t0, s0 = state_at(veh, site, el, az, t_drop, wind=wind)
    body = vehicle.BallisticBody(1.0, beta, "stage")
    tr = trajectory.propagate(body, t0, s0, dt_atm=0.05, dt_vac=0.5, record_every=50, wind=wind)
    return tr.lat[-1], tr.lon[-1], tr.t[-1], tr.apogee()[1]
