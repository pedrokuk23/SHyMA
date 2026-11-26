"""Two-parameter calibration of the VSB-30 model on a single flight.

The effective specific impulse (common to S31 and S30) and the effective
(equivalent constant-thrust) burn time of the S31 booster are adjusted so
that the model reproduces the measured apogee and ground range of one
flight. By default the first VSB-30 flight is used (V01, Alcantara, 392 kg
payload, 84 deg elevation, 45 deg azimuth, apogee 236.0 km, range 189.0 km;
Garcia et al., JATM, 2011). The remaining flights are held-out validation cases
(scripts/run_vsb30_validation.py, which also performs a leave-one-out study
by calibrating on each flight in turn).

The problem is square (two parameters, two observables), so it is solved
with a damped Newton iteration on a finite-difference Jacobian. The Jacobian
is also returned, to propagate the observation uncertainty of the
calibration flight into the parameters.
"""
import json
import os

from . import analysis, earth, flight, vehicle

V01 = {"flight": "V01", "payload_kg": 392.0, "elevation": 84.0, "azimuth": 45.0,
       "site": "CLA", "apogee_km": 236.0, "range_km": 189.0}
MODEL = "5dof"
DT = 0.01


def fly(isp, b1, f=V01, dt=DT):
    v = vehicle.vsb30(f["payload_kg"], isp=isp, s31_burn=b1)
    site = earth.STATIONS[f["site"]]
    tr = flight.fly(v, site, f["elevation"], f["azimuth"], 60.0, record_every=5, model=MODEL, dt=dt)
    return tr.apogee()[1] / 1000.0, analysis.downrange_km(tr, site)[-1]


def jacobian(isp, b1, f=V01, d_isp=1.0, d_b1=0.2):
    a0, r0 = fly(isp, b1, f)
    a1, r1 = fly(isp + d_isp, b1, f)
    a2, r2 = fly(isp, b1 + d_b1, f)
    return (a0, r0), [[(a1 - a0) / d_isp, (a2 - a0) / d_b1], [(r1 - r0) / d_isp, (r2 - r0) / d_b1]]


def calibrate(f=V01, x0=(252.0, 12.0), tol_km=0.02, max_iter=12):
    isp, b1 = x0
    hist = []
    for _ in range(max_iter):
        (a, r), J = jacobian(isp, b1, f)
        ea, er = f["apogee_km"] - a, f["range_km"] - r
        hist.append((isp, b1, a, r))
        if abs(ea) < tol_km and abs(er) < tol_km:
            break
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        dx = (J[1][1] * ea - J[0][1] * er) / det
        dy = (-J[1][0] * ea + J[0][0] * er) / det
        sc = max(abs(dx) / 10.0, abs(dy) / 3.0, 1.0)
        isp, b1 = isp + dx / sc, min(max(b1 + dy / sc, 4.0), 20.0)
    (a, r), J = jacobian(isp, b1, f)
    return {"isp": isp, "b1": b1, "apogee_km": a, "range_km": r, "jacobian": J,
            "iterations": len(hist)}


def parameter_covariance(J, sig_apo_km, sig_rng_km, rho=0.0):
    """Covariance of (Isp, t_b1) induced by the observation covariance of the
    calibration flight (standard deviations of apogee and range and their
    correlation rho, from the Monte Carlo model): Jinv S Jinv^T."""
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    Ji = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]]
    c = rho * sig_apo_km * sig_rng_km
    S = [[sig_apo_km ** 2, c], [c, sig_rng_km ** 2]]
    C = [[sum(Ji[i][k] * S[k][l] * Ji[j][l] for k in range(2) for l in range(2)) for j in range(2)]
         for i in range(2)]
    return C


if __name__ == "__main__":
    c = calibrate()
    out = {"isp_vsb30": round(c["isp"], 3), "s31_burn": round(c["b1"], 4), "model": MODEL,
           "earth_model": earth.MODEL, "calibration_flight": "V01",
           "exit_area_m2": list(vehicle.EXIT_AREA),
           "target_apogee_km": V01["apogee_km"], "target_range_km": V01["range_km"],
           "apogee_check_km": c["apogee_km"], "range_check_km": c["range_km"],
           "jacobian_km_per_s": c["jacobian"], "iterations": c["iterations"]}
    with open(os.path.join(os.path.dirname(__file__), "calibration.json"), "w") as fo:
        json.dump(out, fo, indent=2)
    print(out)
