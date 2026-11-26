"""Code verification.

1) RK4 integrator against the analytic Keplerian solution (drag and J2 off,
   spherical Earth): position error at a fixed time (t = 200 s) for several
   steps, and the observed order of convergence; apogee radius error.
2) Step-size convergence of the nominal VSB-30 point-mass flight and of the
   5-DOF ascent + 3-DOF flight used in the case study (impact point, apogee
   and test-window dwell time). Impact, apogee and window boundaries are
   located by interpolation inside the step.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from shyma import earth, flight, trajectory, vehicle  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
res = {}


def kepler_position(r0, v0, t, mu):
    """Analytic two-body position at time t (elliptic, planar x-z), from
    Lagrange f and g functions in the eccentric-anomaly increment dE."""
    rx, rz = r0
    vx, vz = v0
    r = math.hypot(rx, rz)
    a = 1.0 / (2.0 / r - (vx * vx + vz * vz) / mu)
    n = math.sqrt(mu / a ** 3)
    e_cos = 1 - r / a                        # e cos E0
    e_sin = (rx * vx + rz * vz) / math.sqrt(mu * a)   # e sin E0
    dE = n * t
    for _ in range(60):
        F = dE - e_cos * math.sin(dE) + e_sin * (1 - math.cos(dE)) - n * t
        dE -= F / (1 - e_cos * math.cos(dE) + e_sin * math.sin(dE))
    f = 1 - a / r * (1 - math.cos(dE))
    g = t - (dE - math.sin(dE)) / n
    return f * rx + g * vx, f * rz + g * vz


# 1) Keplerian check
earth.set_model("sphere", j2=False, rotation=False)
body = vehicle.BallisticBody(100.0, 1e15)
r0 = earth.RE + 100e3
v0, g0 = 2500.0, math.radians(60.0)
s0 = [r0, 0.0, 0.0, v0 * math.sin(g0), 0.0, v0 * math.cos(g0)]
eps = 0.5 * v0 ** 2 - earth.MU / r0
hmom = r0 * v0 * math.cos(g0)
a = -earth.MU / (2 * eps)
e = math.sqrt(1 + 2 * eps * hmom ** 2 / earth.MU ** 2)
ra_analytic = a * (1 + e)
T_CHECK = 200.0
xa, za = kepler_position((s0[0], s0[2]), (s0[3], s0[5]), T_CHECK, earth.MU)
errs = []
for dt in (20.0, 10.0, 5.0, 2.5):
    tr = trajectory.propagate(body, 0.0, s0, dt_atm=dt, dt_vac=dt, stop_alt=99e3, t_end=T_CHECK)
    x, y, z = tr.final[1][:3]
    err = math.sqrt((x - xa) ** 2 + y ** 2 + (z - za) ** 2)
    tr2 = trajectory.propagate(body, 0.0, s0, dt_atm=dt, dt_vac=dt, stop_alt=99e3)
    ra_num = tr2.apogee()[1] + earth.RE
    res["kepler_dt_%.2f" % dt] = {"pos_err_m_t200": err, "ra_err_m": ra_num - ra_analytic}
    errs.append((dt, err))
res["kepler_observed_order"] = [math.log(e1 / e2) / math.log(d1 / d2)
                                for (d1, e1), (d2, e2) in zip(errs[:-1], errs[1:])]
res["kepler_ra_analytic_km"] = (ra_analytic - earth.RE) / 1000.0
earth.reset_model()

# 2) VSB-30 point mass (3-DOF) and case-study driver (5-DOF + 3-DOF)
site = earth.STATIONS["CLA"]
for tag, mk in (("vsb30", "3dof"), ("fly5", "5dof")):
    out = {}
    for dt in (0.1, 0.05, 0.02, 0.01, 0.005):
        v = vehicle.vsb30(400.0, t_payload_sep=60.0)
        tr = flight.fly(v, site, 84.0, 45.0, 60.0, dt=dt, record_every=10, model=mk)
        out[dt] = (tr.lat[-1], tr.lon[-1], tr.apogee()[1], tr.t[-1])
    ref = out[0.005]
    for dt in (0.1, 0.05, 0.02, 0.01):
        o = out[dt]
        res["%s_dt_%.3f" % (tag, dt)] = {
            "impact_diff_m": earth.great_circle_km(o[0], o[1], ref[0], ref[1]) * 1000.0,
            "apogee_diff_m": o[2] - ref[2], "t_impact_diff_s": o[3] - ref[3]}
    res["%s_ref" % tag] = {"apogee_km": ref[2] / 1e3, "t_impact": ref[3]}

# 3) window dwell time convergence (case-study concepts)
import run_casestudy as cs  # noqa: E402
for c in "ABCD":
    d = {}
    for dt in (0.05, 0.02, 0.01, 0.005):
        v, tr = cs.fly(c, 400.0, cs.NOMINAL[c], dt=dt, rec=1)
        d[dt] = cs.windows(c, v, tr)[0]["duration"]
    res["window_%s" % c] = {"%.3f" % k: x for k, x in d.items()}
with open(os.path.join(OUT, "verification.json"), "w") as f:
    json.dump(res, f, indent=2)
print(json.dumps(res, indent=2))
