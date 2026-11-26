"""Concept D: recovery of the stack from an attitude error at S30 ignition.

The S30 of Concept D ignites after a long unpowered coast at low dynamic
pressure. An attitude error at ignition (angle 0, 2, 5 deg, two phases) is
applied to the 5-DOF stack, with the generic restoring derivatives of the
paper and with C_N_alpha halved (weaker static stability), and the
angle-of-attack envelope (time after which it stays below 0.5 deg), the
largest q*alpha of the S30 burn, the burnout state and the dwell are recorded.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_casestudy as cs  # noqa: E402
import run_montecarlo as mc  # noqa: E402
from shyma import flight  # noqa: E402

ANGLES = (0.0, 2.0, 5.0)
PHASES = (0.0, math.pi / 2)
SCALES = (1.0, 0.5)


def case(angle, phase, scale):
    base = dict(flight.VSB30_AERO)
    flight.VSB30_AERO = dict(base, cn_alpha=base["cn_alpha"] * scale)
    try:
        v, tr = cs.fly("D", 400.0, cs.NOMINAL["D"], rec=1,
                       ign_att=None if angle == 0.0 else (angle, phase))
    finally:
        flight.VSB30_AERO = base
    t_ign = v.stages[1].t_ign
    w = cs.windows("D", v, tr)[0]
    al = tr.alpha
    ks = [k for k in range(min(len(al), len(tr.t))) if t_ign <= tr.t[k] <= v.t_burnout]
    amax = max(al[k] for k in ks) if ks else 0.0
    # envelope criterion: time after ignition beyond which |alpha| stays below 0.5 deg until burnout
    t_half = 0.0
    for k in ks:
        if al[k] > 0.5:
            t_half = tr.t[k] - t_ign
    qa = max(tr.q[k] * al[k] for k in ks) / 1e3 if ks else 0.0
    # angle-of-attack amplitude 2 s and 5 s after ignition (max over the following 1.5 s, about one period)
    def amp(t0):
        w = [al[k] for k in ks if t_ign + t0 <= tr.t[k] <= t_ign + t0 + 1.5]
        return max(w) if w else None
    s = cs.summarize("D", v, tr, light=True)
    return {"angle_deg": angle, "phase": phase, "cn_scale": scale, "t_ign": t_ign,
            "alpha_max_deg": amax, "t_alpha_below_0p5_s": t_half, "qalpha_max_kpa_deg": qa,
            "alpha_amp_2s_deg": amp(2.0), "alpha_amp_5s_deg": amp(5.0),
            "alpha_at_burnout_deg": al[ks[-1]] if ks else None,
            "window_s": w["duration"], "burnout_mach": s["burnout_mach"],
            "burnout_alt_km": s["burnout_alt_km"],
            "tau_use": cs.useful_window(tr, w, ant_loss_db=mc.ANT_LOSS_DB, fade_db=mc.FADE_DB,
                                        mask_deg=mc.ELEV_MASK_DEG)}


if __name__ == "__main__":
    res = []
    for sc in SCALES:
        for a in ANGLES:
            for ph in (PHASES if a > 0 else (0.0,)):
                r = case(a, ph, sc)
                print(r, flush=True)
                res.append(r)
    json.dump({"cases": res}, open(os.path.join(os.path.dirname(__file__), "..", "results", "ignition_check.json"), "w"),
              indent=1)
