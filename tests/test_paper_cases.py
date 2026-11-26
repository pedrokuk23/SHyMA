"""pytest: every quick paper case (cases/*.json) reproduces the published values."""
import glob
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "cases"))
import run_case  # noqa: E402

CASES = sorted(glob.glob(os.path.join(ROOT, "cases", "*.json")))
QUICK = [f for f in CASES if not json.load(open(f)).get("slow")]


@pytest.mark.parametrize("path", QUICK, ids=[os.path.basename(f)[:-5] for f in QUICK])
def test_case(path):
    ok, rows, _ = run_case.run(json.load(open(path)))
    assert ok, [r for r in rows if not r[3]]
