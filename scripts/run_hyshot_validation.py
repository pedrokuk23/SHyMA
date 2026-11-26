"""Validation against the HyShot 2 flight (30 July 2002, Woomera).

Re-entry reconstruction from the state published by Hass, Smart & Paull
(AIAA 2005-3354, NASA NTRS 20050215040, Trajectory Re-construction section):
the re-entry phase starts at t = 508 s, 99.25 km altitude, 2172.9 m/s and a
-69.2 deg flight-path angle (obtained there by integrating the axial
accelerometer and matching the Pitot pressure history). The same paper gives,
for the end of its space phase, t = 512.8 s at the same altitude with
2185.6 m/s and -69.6 deg; that state is run as a sensitivity case. Heading and
position are not published (a 330 deg heading over Woomera is used; both have
a negligible effect), and the velocity is taken as Earth-relative.
Test window from Smart, Hass & Paull (AIAA J., 2006): ~3 s between 35 and
29 km altitude, Mach 7.75-8.0, dynamic pressure 21.0-58.5 kPa; window opens
at t = 537.5 s. The re-entry body ballistic coefficient is not published and
is treated as an uncertain parameter (sweep).
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from shyma import analysis, earth, trajectory, vehicle  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results")

ENTRY = {"t": 508.0, "alt_km": 99.25, "v": 2172.9, "fpa": -69.2}       # re-entry initial state
SPACE_END = {"t": 512.8, "alt_km": 99.25, "v": 2185.6, "fpa": -69.6}   # end of the space phase
APOGEE_KM = 328.27
LAT, LON = earth.STATIONS["Woomera"][:2]
HEADING = 330.0

FLIGHT = {"alt_km": [35.0, 29.0], "mach": [7.75, 8.0], "q_kpa": [21.0, 58.5],
          "duration_s": 3.0, "t_window_open": 537.5}


def run(beta=5000.0, rho_factor=1.0, heading=HEADING, atm=None, st=ENTRY):
    body = vehicle.BallisticBody(110.0, beta, "HyShot2")
    g, hd = math.radians(st["fpa"]), math.radians(heading)
    vh = st["v"] * math.cos(g)
    s0 = trajectory.state_from_ned(LAT, LON, st["alt_km"] * 1000.0, vh * math.cos(hd), vh * math.sin(hd),
                                   -st["v"] * math.sin(g), st["t"])
    return trajectory.propagate(body, st["t"], s0, rho_factor=rho_factor, dt_atm=0.01,
                                dt_vac=0.1, stop_alt=15e3, atm=atm)


def window_summary(tr, h_top=35.0, h_bot=29.0):
    a35 = analysis.state_at_altitude(tr, h_top)
    a29 = analysis.state_at_altitude(tr, h_bot)
    return {"t_35km": a35["t"], "t_29km": a29["t"], "duration_35_29": a29["t"] - a35["t"],
            "mach_35": a35["mach"], "mach_29": a29["mach"],
            "q_35_kpa": a35["q"] / 1e3, "q_29_kpa": a29["q"] / 1e3,
            "v_35": a35["v"], "fpa_35": a35["fpa"]}


if __name__ == "__main__":
    res = {"flight": FLIGHT, "entry_state": ENTRY, "space_end_state": SPACE_END}
    res["nominal"] = window_summary(run())
    # the 2005 report associates the measured q rise with 34.9 -> 28.1 km (the 2006 paper rounds to 35 -> 29 km)
    res["ref_alt_2005"] = window_summary(run(), 34.9, 28.1)
    # sensitivity cases
    sens = {}
    for b in (1500.0, 5000.0, 15000.0):
        sens["beta_%d" % b] = window_summary(run(beta=b))
    for f in (0.85, 1.15):
        sens["rho_x%.2f" % f] = window_summary(run(rho_factor=f))
    sens["space_end_state"] = window_summary(run(st=SPACE_END))
    # start from the reconstructed apogee (configuration of an earlier version of the analysis): horizontal
    # speed equal to that of the published re-entry state (771.6 m/s); the flight time of the apogee is
    # not published, so only the Mach number and q at 35 and 29 km are meaningful, not the timing
    vh = ENTRY["v"] * math.cos(math.radians(ENTRY["fpa"]))
    sens["apogee_start"] = window_summary(run(st={"t": 278.0, "alt_km": APOGEE_KM, "v": vh, "fpa": 0.0}))
    sens["heading_150"] = window_summary(run(heading=150.0))
    # site- and date-specific atmosphere (NRLMSISE-00, Woomera, 30 July 2002;
    # F10.7 = 180 sfu and Ap = 10 assumed for the 2002 solar maximum)
    from shyma.atmosphere import msis_atmosphere
    msis = msis_atmosphere(LAT, LON, "2002-07-30T12:00", f107=180.0, f107a=180.0, ap=10.0)
    sens["nrlmsise00"] = window_summary(run(atm=msis))
    # static pressure implied by the flight data, p = 2 q / (gamma M^2), vs US76
    from shyma.atmosphere import atmosphere as us76
    res["implied_pressure"] = {
        "35km": {"flight_Pa": 2 * 21.0e3 / (1.4 * 7.75 ** 2), "us76_Pa": us76(35.0)[1],
                 "msis_Pa": msis(35.0)[1]},
        "29km": {"flight_Pa": 2 * 58.5e3 / (1.4 * 8.0 ** 2), "us76_Pa": us76(29.0)[1],
                 "msis_Pa": msis(29.0)[1]}}
    res["sensitivity"] = sens
    # relative errors (nominal) against flight
    n = res["nominal"]
    res["errors"] = {
        "mach_35_pct": 100 * (n["mach_35"] - 7.75) / 7.75,
        "mach_29_pct": 100 * (n["mach_29"] - 8.0) / 8.0,
        "q_35_pct": 100 * (n["q_35_kpa"] - 21.0) / 21.0,
        "q_29_pct": 100 * (n["q_29_kpa"] - 58.5) / 58.5,
        "duration_s": n["duration_35_29"] - 3.0,
        "t_open_s": n["t_35km"] - 537.5,
    }
    # store the nominal descent for plotting
    tr = run()
    res["profile"] = {"t": tr.t[::20], "alt_km": [a / 1e3 for a in tr.alt[::20]],
                      "mach": tr.mach[::20], "q_kpa": [q / 1e3 for q in tr.q[::20]]}
    with open(os.path.join(OUT, "hyshot_validation.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps({k: v for k, v in res.items() if k != "profile"}, indent=2))
