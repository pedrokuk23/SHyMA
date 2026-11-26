"""Geometry-based stability derivatives: exact limits and the inverse dynamic check."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import run_stability_check as chk  # noqa: E402
from shyma import flight  # noqa: E402


def test_exact_limits():
    r = chk.limits()
    assert abs(r["cone_cna"] - 2.0) < 1e-12 and abs(r["cone_xcp_over_L"] - 2.0 / 3.0) < 1e-12
    assert abs(r["diederich_large_A_over_2pi"] - 1.0) < 2e-3
    assert abs(r["diederich_small_A_over_piA2"] - 1.0) < 1e-5
    assert abs(r["supersonic_large_A_times_beta_over_4"] - 1.0) < 1e-3
    assert abs(r["sonic_le_jump"]) < 1e-4
    assert abs(r["kwb_0"] - 1.0) < 1e-6 and abs(r["kwb_1"] - 2.0) < 5e-3
    assert r["kwb_plus_kbw_minus_1ptau2"] < 1e-12
    assert abs(r["E_05"] - r["E_05_ref"]) < 1e-12


def test_inverse_dynamic_check():
    r = chk.inverse_check()
    assert abs(r["err_cm_alpha_pct"]) < 0.1 and abs(r["err_cmq_pct"]) < 0.5


def test_default_is_generic():
    assert flight.AERO_FACTORY is None


def test_drag_building_blocks():
    r = chk.drag_checks()
    assert abs(r["pm_nu_M2_deg"] - r["pm_nu_M2_ref"]) < 0.01
    assert abs(r["oblique_M2_d10_beta_deg"] - r["oblique_beta_ref"]) < 0.02
    assert abs(r["oblique_M2_d10_p21"] - r["oblique_p21_ref"]) < 1e-3
    assert abs(r["cone10_M2_shock_deg"] - r["cone10_shock_ref"]) < 0.1
    assert 0.94 < r["cone10_M6_cp_over_hsd"] < 1.03 and 0.94 < r["cone5_M8_cp_over_hsd"] < 1.03
    assert abs(r["cf_KS_Re1e7"] - r["cf_KS_ref"]) < 2e-6
    assert 0.75 < r["cf_VD_ratio_M2_Re1e7"] < 0.88
