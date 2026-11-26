"""Regression tests for three code defects found in review (they do not affect the paper's cases)."""
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from shyma import earth, trajectory5 as t5, vehicle  # noqa: E402
from shyma.trajectory import Env  # noqa: E402
import launch_optimizer as lo  # noqa: E402


def test_dropped_stage_gives_no_thrust():
    v = vehicle.vsb30(400.0)
    s30 = [s for s in v.stages if s.name == "S30"][0]
    s30.t_drop = s30.t_ign + 0.5 * s30.burn_time          # drop the stage in the middle of its burn
    t = s30.t_ign + 0.75 * s30.burn_time
    assert v.thrust(t) == 0.0


def test_feasible_setting_ranks_above_infeasible():
    base = {"threshold_met": True, "closeness": 0.0, "q_max_ascent_kpa": 100.0}
    feas = dict(base, threshold=4.0, impact_range_km=299.9)
    infeas = dict(base, threshold=4.1, impact_range_km=300.001)
    cons = {"range_max": 300.0}
    sf, ok_f = lo.score(feas, "max_dwell", "threshold", cons)
    si, ok_i = lo.score(infeas, "max_dwell", "threshold", cons)
    assert ok_f and not ok_i and sf > si


def test_5dof_recorded_acceleration_matches_dynamics():
    veh = vehicle.vsb30(400.0)
    aero = t5.Aero5(0.557)
    env = Env()
    r0 = earth.RE + 20e3
    a0 = math.radians(3.0)
    s = [r0, 0.0, 0.0, 0.0, 900.0, 0.0, 0.0, math.cos(a0), math.sin(a0), 0.0, 0.0, 0.0]
    tr = t5.propagate5(veh, aero, 0.0, s, t_end=0.02, dt_atm=0.01, record_every=1)
    d5 = t5._deriv5(0.0, s, veh, aero, env)
    g = earth.gravity(*s[:3])
    acc = math.sqrt(sum((d5[3 + i] - g[i]) ** 2 for i in range(3)))
    assert abs(tr.acc[0] - acc) < 1e-9 * max(1.0, acc)
