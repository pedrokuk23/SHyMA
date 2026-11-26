"""Useful test time and three-state requirement report."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

import run_casestudy as cs  # noqa: E402
import launch_optimizer as lo  # noqa: E402


def _flight(c, el):
    v, tr = cs.fly(c, 400.0, el, dt=0.02, rec=2)
    return v, tr, cs.windows(c, v, tr)[0]


def test_useful_time_bounded_by_window():
    v, tr, w = _flight("B", 84.0)
    rec = cs.useful_window(tr, w, need_link=False)
    live = cs.useful_window(tr, w)
    assert abs(rec - w["duration"]) < 1e-9          # recording: the whole flight window
    assert 0.0 <= live <= rec + 1e-9                 # live data: a continuous part of it
    assert cs.useful_window(tr, w, t_ready=1.0, need_link=False) == w["duration"] - 1.0


def test_low_elevation_reentry_loses_cla_link():
    v, tr, w = _flight("B", 66.0)
    assert w["met"] and w["duration"] > 3.0
    assert cs.useful_window(tr, w) == 0.0            # window far downrange, below the CLA mask


def test_requirement_report_has_three_states():
    v, tr, w = _flight("B", 84.0)
    wo = lo.window_objects("B", v, cs.DSEP)["threshold"]
    d = lo.diagnose(wo, tr, "B")
    states = {i["status"] for i in d["items"]}
    assert states == {"met", "not evaluated"}
    assert all(i["ok"] is None for i in d["items"] if i["status"] == "not evaluated")
    assert "status" in lo.diagnose(wo, tr)["items"][0]
