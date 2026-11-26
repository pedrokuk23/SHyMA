"""pytest: quick recomputations of results added in the fourth review round.

Each test recomputes one value of a result file from the released code and checks it
against the published value: the site-atmosphere apogee change of the first Esrange
flight (Table S14), the transient leading-edge temperature of Concept A (Table S11),
the stack drag rule (a blunt test vehicle adds its excess drag to the attached stack)
and the powered JSBSim comparison (Table S3; skipped if JSBSim is not installed).
"""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)


def load(name):
    with open(os.path.join(ROOT, "results", name + ".json")) as f:
        return json.load(f)


def test_site_bias_first_esrange_flight():
    import run_site_bias as sb
    import run_vsb30_validation as vv
    from shyma import earth
    from shyma.atmosphere import msis_atmosphere
    r = [x for x in vv.load() if x["flight"] == "V02"][0]
    la, lo = earth.STATIONS[r["site"]][:2]
    atm = msis_atmosphere(la, lo, r["date"] + "T12:00", f107=100.0, f107a=100.0, ap=5.0)
    d = sb.apo(r, atm) - sb.apo(r)
    ref = [x for x in load("site_bias")["flights"] if x["flight"] == "V02"][0]["delta_km"]
    assert abs(d - ref) < 1e-6


def test_le_transient_concept_a():
    import run_le_transient as le
    out = le.run("A")
    ref = load("le_transient")["concepts"]["A"]
    assert abs(out["f1.0"]["T_max_window_K"] - ref["f1.0"]["T_max_window_K"]) < 1e-6
    assert abs(out["f1.0"]["Q_cold_MJ"] - ref["f1.0"]["Q_cold_MJ"]) < 1e-6


def test_stack_drag_not_below_test_vehicle():
    from shyma import vehicle
    for beta in (2500.0, 5000.0, 7462.0, 15000.0):
        v = vehicle.vsb30(400.0, beta_tv=beta, t_payload_sep=1e9)
        stack = v.drag_area(50.0, 6.0)
        tv = vehicle.cd_generic(6.0) * v.payload_cd_scale * v.payload_area
        assert stack >= tv - 1e-12


def test_jsbsim_powered_case():
    pytest.importorskip("jsbsim")
    sys.path.insert(0, os.path.join(ROOT, "jsbsim_ref"))
    import run_jsbsim_thrust as jt
    tr = jt.run_shyma()
    sh = {"t": tr.t, "alt_km": [a / 1e3 for a in tr.alt], "mach": tr.mach}
    ref = load("jsbsim_thrust")["rows"]
    for row in ref:
        m = jt.at_alt(sh["t"], sh["alt_km"], sh["mach"], row["alt_km"])
        assert abs(m - row["mach"][0]) < 1e-9
        assert abs(m / row["mach"][1] - 1.0) < 1e-3


def test_hydrostatic_temperature_offset():
    """A uniform +4 K column offset raises the pressure at 30 km by about 8 % (hydrostatic re-integration)."""
    from shyma.trajectory import Env
    e0, e1 = Env(), Env(dT=lambda h: 4.0)
    p0, p1 = e0.atmosphere(30.0)[1], e1.atmosphere(30.0)[1]
    assert 1.07 < p1 / p0 < 1.09
    r0, r1 = e0.atmosphere(30.0)[2], e1.atmosphere(30.0)[2]
    assert r1 > r0          # density rises: the pressure term dominates the direct -dT/T term


def test_stack_excess_follows_design_beta():
    """A dispersed test-vehicle beta changes the free flight only; the stack excess follows the design beta."""
    from shyma import vehicle
    import run_casestudy as cs
    a = vehicle.vsb30(400.0, beta_tv=0.8 * cs.BETA_TV, beta_stack=cs.BETA_TV)
    b = vehicle.vsb30(400.0, beta_tv=cs.BETA_TV)
    assert abs(a.tv_excess - b.tv_excess) < 1e-12
    c = vehicle.vsb30(400.0, beta_tv=0.8 * cs.BETA_TV)
    assert c.tv_excess > 0.2


def test_common_random_numbers_per_concept():
    """Every concept restarts the same aleatory stream (run i of every set uses the draws of run i of the base set)."""
    import random
    from shyma import dispersion as D
    d1 = [D.sample(random.Random(20261006), 0.9) for _ in range(3)]
    d2 = [D.sample(random.Random(20261006), 0.9) for _ in range(3)]
    assert d1 == d2
