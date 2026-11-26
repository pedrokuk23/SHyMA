"""Monte Carlo dispersion models shared by the wind calibration, the VSB-30
validation and the case-study Monte Carlo.

1-sigma values (Gaussian). No vehicle-specific dispersion data are public
for the VSB-30, so they are engineering estimates of the order used in
conceptual dispersion analyses of solid sounding rockets; Section VI-F of
the paper reports the sensitivity of the outputs to doubling all of them.

  imp   stage total impulse (propellant mass at fixed Isp)        1 %
  burn  stage burn time at constant impulse (propellant temp.)    2 %
  dry   stage dry mass                                            1 %
  pl    payload mass                                              0.5 %
  cd    axial-drag scale                                          5 %
  rho   atmospheric pressure/density scale (whole profile)        5 %
  dT    atmospheric temperature offset, altitude-correlated       4 K
  el    launcher elevation (rail, tip-off)                        0.2 deg
  az    launcher azimuth                                          0.3 deg
  beta  test-vehicle ballistic coefficient (case study only)      10 %
  dsep  separation/shroud/settling delay (Concept A only)         0.3 s
  coast interstage coast / S30 ignition timer (Concept D only)   0.1 s
  wind  residual-wind profile, amplitude k_w (lumped residual
        dispersion parameter, calibrated: run_wind_calibration.py)

Thrust and fin misalignment, dynamic imbalance and roll coupling are not
modelled explicitly (no roll dynamics); their contribution is lumped into
k_w, which is therefore specific to the calibration configuration.
"""
import math
import random

from . import analysis, earth, flight, vehicle, wind

SIG = {"imp": 0.01, "burn": 0.02, "dry": 0.01, "pl": 0.005, "cd": 0.05, "rho": 0.05,
       "dT": 4.0, "el": 0.2, "az": 0.3, "beta": 0.10, "dsep": 0.3, "coast": 0.1}


def sample(rng, k_wind, scale=1.0):
    s = {k: v * scale for k, v in SIG.items()}
    return {"imp1": rng.gauss(1, s["imp"]), "imp2": rng.gauss(1, s["imp"]),
            "burn1": rng.gauss(1, s["burn"]), "burn2": rng.gauss(1, s["burn"]),
            "dry": rng.gauss(1, s["dry"]), "pl": rng.gauss(1, s["pl"]),
            "cd": rng.gauss(1, s["cd"]), "rho": rng.gauss(1, s["rho"]),
            "dT0": rng.gauss(0, s["dT"]), "dT1": rng.gauss(0, s["dT"]),
            "del": rng.gauss(0, s["el"]), "daz": rng.gauss(0, s["az"]),
            "beta": rng.gauss(1, s["beta"]), "dsep": rng.gauss(0, s["dsep"]),
            "dcoast": rng.gauss(0, s["coast"]),
            "wseed": rng.randrange(1 << 30), "k_wind": k_wind}


def wind_of(d):
    if d["k_wind"] <= 0.0:
        return None
    return wind.residual_wind(random.Random(d["wseed"]), d["k_wind"])


class TempOffset:
    """Altitude-correlated temperature offset: a lower (0-15 km) and an upper
    (above 30 km) random value, blended linearly between 15 and 30 km."""

    def __init__(self, dt0, dt1):
        self.a, self.b = dt0, dt1

    def __call__(self, alt_m):
        c = min(max((alt_m - 15e3) / 15e3, 0.0), 1.0)
        return ((1 - c) * self.a + c * self.b) / math.sqrt((1 - c) ** 2 + c ** 2)


def temp_of(d):
    return TempOffset(d["dT0"], d["dT1"])


def _xy(lat0, lon0, lats, lons):
    xs = [(lo - lon0) * math.cos(math.radians(lat0)) * earth.WGS84_A * math.pi / 180 / 1e3 for lo in lons]
    ys = [(la - lat0) * earth.WGS84_A * math.pi / 180 / 1e3 for la in lats]
    return xs, ys


_GAUSS_U = None


def _gauss_radius(l1, l2, p):
    """Radius containing probability p of a centred bivariate normal with
    variances l1, l2 (fixed pseudo-random sample with a fixed seed, deterministic)."""
    global _GAUSS_U
    if _GAUSS_U is None:
        r = random.Random(12345)
        _GAUSS_U = [(r.gauss(0, 1), r.gauss(0, 1)) for _ in range(200000)]
    a, b = math.sqrt(l1), math.sqrt(max(l2, 0.0))
    rr = sorted(math.hypot(a * u, b * v) for u, v in _GAUSS_U)
    return rr[int(p * len(rr)) - 1]


def ellipse(lat0, lon0, lats, lons):
    xs, ys = _xy(lat0, lon0, lats, lons)
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs) / (n - 1)
    syy = sum((y - my) ** 2 for y in ys) / (n - 1)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (n - 1)
    tr, det = sxx + syy, sxx * syy - sxy ** 2
    l1 = tr / 2 + math.sqrt(max(tr ** 2 / 4 - det, 0.0))
    l2 = tr / 2 - math.sqrt(max(tr ** 2 / 4 - det, 0.0))
    ang = 0.5 * math.degrees(math.atan2(2 * sxy, sxx - syy))
    r = sorted(math.hypot(x - mx, y - my) for x, y in zip(xs, ys))
    p3 = 1.0 - math.exp(-4.5)
    return {"sigma_major_km": math.sqrt(l1), "sigma_minor_km": math.sqrt(max(l2, 0.0)),
            "angle_from_east_deg": ang, "r99_km": r[int(0.99 * n) - 1],
            "r3sig_empirical_km": r[int(p3 * n) - 1],
            "r3sig_gauss_km": _gauss_radius(l1, l2, p3),
            "r_max_km": r[-1], "x_km": xs, "y_km": ys, "mean_xy": (mx, my)}


def clopper_pearson(k, n, conf=0.95):
    from scipy.stats import beta
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return float(lo), float(hi)


def percentile(xs, p):
    s = sorted(xs)
    if not s:
        return float("nan")
    i = p * (len(s) - 1)
    k = int(i)
    return s[k] if k + 1 >= len(s) else s[k] + (i - k) * (s[k + 1] - s[k])


def bootstrap_ci(xs, stat, n_boot=1000, conf=0.95, seed=1):
    r = random.Random(seed)
    n = len(xs)
    vals = sorted(stat([xs[r.randrange(n)] for _ in range(n)]) for _ in range(n_boot))
    a = (1 - conf) / 2
    return vals[int(a * n_boot)], vals[int((1 - a) * n_boot) - 1]



def vsb30_dispersed(payload, d, **kw):
    return vehicle.vsb30(payload * d["pl"], cd_scale=d["cd"], payload_cd_scale=d["cd"],
                         impulse_factor=(d["imp1"], d["imp2"]), dry_factor=d["dry"],
                         burn_factor=(d["burn1"], d["burn2"]), **kw)


def one_flight(args):
    """Dispersed VSB-30 flight for a configuration (payload, el, az, site)."""
    (payload, el, az, site_name), d = args
    v = vsb30_dispersed(payload, d)
    site = earth.STATIONS[site_name]
    kick = None if "tip" not in d else (None, d["tip"][0], d["tip"][1])
    tr = flight.fly(v, site, el + d["del"], az + d["daz"], 60.0, wind=wind_of(d),
                    rho_factor=d["rho"], dT=temp_of(d), dt=0.02, record_every=25, att_kick=kick)
    return {"impact_lat": tr.lat[-1], "impact_lon": tr.lon[-1],
            "range_km": analysis.downrange_km(tr, site)[-1], "apogee_km": tr.apogee()[1] / 1e3}
