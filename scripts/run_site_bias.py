"""Launch-site atmosphere and the CLA/Esrange split of the held-out errors.

The booster is calibrated on V01 (CLA, tropical) with the US Standard
Atmosphere 1976 and predicts six Esrange (subarctic) flights. Here every
flight is also flown with NRLMSISE-00 for its site and launch date (12 UTC,
F10.7 = F10.7a = 100, Ap = 5). Because a recalibration on V01 with
NRLMSISE-00 absorbs the V01 change, the first-order site-atmosphere effect on
each held-out prediction is Delta_i - Delta_V01, where Delta is the apogee
change from US76 to NRLMSISE-00 with the same booster parameters.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import run_vsb30_validation as vv  # noqa: E402
from shyma import earth, flight, vehicle  # noqa: E402
from shyma.atmosphere import msis_atmosphere  # noqa: E402


def apo(r, atm=None):
    v = vehicle.vsb30(float(r["payload_kg"]))
    tr = flight.fly(v, earth.STATIONS[r["site"]], float(r["elevation_deg"]), float(r["azimuth_deg"]), 60.0,
                    record_every=5, atm=atm)
    return tr.apogee()[1] / 1e3


if __name__ == "__main__":
    rows = vv.load()
    out = []
    for r in rows:
        la, lo = earth.STATIONS[r["site"]][:2]
        atm = msis_atmosphere(la, lo, r["date"] + "T12:00", f107=100.0, f107a=100.0, ap=5.0)
        a0, a1 = apo(r), apo(r, atm)
        out.append({"flight": r["flight"], "site": r["site"], "actual_km": float(r["apogee_actual_km"]),
                    "us76_km": a0, "msis_km": a1, "delta_km": a1 - a0})
        print(out[-1], flush=True)
    d01 = [x for x in out if x["flight"] == "V01"][0]["delta_km"]
    for x in out:
        x["err_us76_pct"] = 100 * (x["us76_km"] - x["actual_km"]) / x["actual_km"]
        x["err_msis_pct"] = 100 * (x["us76_km"] + x["delta_km"] - d01 - x["actual_km"]) / x["actual_km"]
    ho = [x for x in out if x["flight"] != "V01"]
    res = {"flights": out, "delta_v01_km": d01}
    for site in ("CLA", "Esrange"):
        s = [x for x in ho if x["site"] == site]
        res[site] = {"n": len(s), "mean_err_us76_pct": sum(x["err_us76_pct"] for x in s) / len(s),
                     "mean_err_msis_pct": sum(x["err_msis_pct"] for x in s) / len(s)}
    rms = lambda k: (sum(x[k] ** 2 for x in ho) / len(ho)) ** 0.5  # noqa: E731
    res["rms_us76_pct"], res["rms_msis_pct"] = rms("err_us76_pct"), rms("err_msis_pct")
    # normalised residuals with the site atmosphere (two-sided chi-square against N(0,1))
    from scipy.stats import chi2
    vj = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "vsb30_validation.json")))
    sg = {f["flight"]: f["mc_sigma_apogee_km"] for f in vj["flights"]}
    z = [(x["us76_km"] + x["delta_km"] - d01 - x["actual_km"]) / sg[x["flight"]] for x in ho]
    q = sum(v * v for v in z)
    res["chi2_msis"] = q
    res["p_upper_msis"], res["p_lower_msis"] = float(chi2.sf(q, len(z))), float(chi2.cdf(q, len(z)))
    print({k: v for k, v in res.items() if k != "flights"})
    json.dump(res, open(os.path.join(os.path.dirname(__file__), "..", "results", "site_bias.json"), "w"), indent=1)
