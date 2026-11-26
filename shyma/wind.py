"""Horizontal wind profiles, residual-wind dispersion model and wind weighting.

A wind profile is a callable alt_m -> (east, north) in m/s (wind *towards*
east/north, i.e. the air velocity). Profiles enter the equations of motion
through the air-relative velocity, so a weathercock-stable unguided rocket
turns into the wind during the low-speed boost.

Residual-wind model (Monte Carlo): after wind weighting, the launcher is set
from a wind measured shortly before lift-off, and the dispersion is driven by
the difference between that measurement and the wind actually flown. This
residual is modelled as two random horizontal vectors (a low-level and an
upper-level component), blended with altitude and scaled by a standard
deviation envelope sigma(h) whose amplitude k is the single parameter of the
model (calibrated against the flight-proven VSB-30 dispersion).
"""
import bisect
import math

from . import earth, flight


class TabulatedWind:
    """Piecewise-linear wind profile; zero above the last point."""

    def __init__(self, alts_m, east, north):
        self.h, self.e, self.n = list(alts_m), list(east), list(north)

    def __call__(self, alt_m):
        h = self.h
        if alt_m <= h[0]:
            return self.e[0], self.n[0]
        if alt_m >= h[-1]:
            return 0.0, 0.0
        k = bisect.bisect_right(h, alt_m)
        w = (alt_m - h[k - 1]) / (h[k] - h[k - 1])
        return (self.e[k - 1] + w * (self.e[k] - self.e[k - 1]),
                self.n[k - 1] + w * (self.n[k] - self.n[k - 1]))


def uniform_wind(speed, from_deg, top_m=20e3):
    """Constant wind of given speed blowing FROM the given azimuth."""
    to = math.radians(from_deg + 180.0)
    return TabulatedWind([0.0, top_m, top_m + 1.0], [speed * math.sin(to)] * 2 + [0.0],
                         [speed * math.cos(to)] * 2 + [0.0])


# Shape of the residual-wind standard deviation (relative units) versus
# altitude: growth through the boundary layer and troposphere, maximum near
# the jet level, decay in the stratosphere, zero above 40 km.
_SIG_H = [0.0, 2e3, 11e3, 20e3, 30e3, 40e3]
_SIG_S = [1.0, 1.5, 3.0, 3.0, 1.0, 0.0]


def sigma_shape(alt_m):
    if alt_m >= _SIG_H[-1]:
        return 0.0
    k = bisect.bisect_right(_SIG_H, alt_m)
    w = (alt_m - _SIG_H[k - 1]) / (_SIG_H[k] - _SIG_H[k - 1])
    return _SIG_S[k - 1] + w * (_SIG_S[k] - _SIG_S[k - 1])


def residual_wind(rng, k_sigma, blend_h=8e3):
    """Random residual wind profile with 1-sigma amplitude k_sigma*sigma_shape(h)."""
    u0 = (rng.gauss(0, 1), rng.gauss(0, 1))
    u1 = (rng.gauss(0, 1), rng.gauss(0, 1))
    hs = [0.0, 1e3, 2e3, 4e3, 6e3, 8e3, 11e3, 15e3, 20e3, 25e3, 30e3, 35e3, 40e3]
    east, north = [], []
    for h in hs:
        c = min(h / blend_h, 1.0)
        s = k_sigma * sigma_shape(h)
        east.append(s * ((1 - c) * u0[0] + c * u1[0]) / math.sqrt((1 - c) ** 2 + c ** 2))
        north.append(s * ((1 - c) * u0[1] + c * u1[1]) / math.sqrt((1 - c) ** 2 + c ** 2))
    return TabulatedWind(hs, east, north)


def impact_point(make_vehicle, site, el, az, wind=None, dt=0.02, t_sep=60.0, **kw):
    v = make_vehicle()
    tr = flight.fly(v, site, el, az, t_sep, wind=wind, dt=dt, record_every=25, **kw)
    return tr.lat[-1], tr.lon[-1], tr


def _en(lat0, lon0, lat, lon):
    x = math.radians(lon - lon0) * math.cos(math.radians(lat0)) * earth.RE
    y = math.radians(lat - lat0) * earth.RE
    return x, y


def wind_weighting(make_vehicle, site, el0, az0, wind, target=None, tol_m=50.0, max_iter=8):
    """Launcher elevation/azimuth that bring the impact point back to the
    no-wind nominal impact (or to `target` = (lat, lon)) in the given wind.
    Newton iteration with a finite-difference Jacobian."""
    if target is None:
        lat_t, lon_t, _ = impact_point(make_vehicle, site, el0, az0, None)
    else:
        lat_t, lon_t = target
    el, az = el0, az0
    hist = []
    for _ in range(max_iter):
        la, lo, _ = impact_point(make_vehicle, site, el, az, wind)
        ex, ey = _en(lat_t, lon_t, la, lo)
        miss = math.hypot(ex, ey)
        hist.append((el, az, miss))
        if miss < tol_m:
            break
        d = 0.2
        la1, lo1, _ = impact_point(make_vehicle, site, el + d, az, wind)
        la2, lo2, _ = impact_point(make_vehicle, site, el, az + d, wind)
        x1, y1 = _en(lat_t, lon_t, la1, lo1)
        x2, y2 = _en(lat_t, lon_t, la2, lo2)
        J = [[(x1 - ex) / d, (x2 - ex) / d], [(y1 - ey) / d, (y2 - ey) / d]]
        det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
        dl = (-J[1][1] * ex + J[0][1] * ey) / det
        da = (J[1][0] * ex - J[0][0] * ey) / det
        # limit the step (the impact is very non-linear near vertical launch)
        sc = max(abs(dl) / 2.0, abs(da) / 10.0, 1.0)
        el = min(max(el + dl / sc, 45.0), 89.5)
        az = (az + da / sc) % 360.0
    return {"elevation": el, "azimuth": az, "miss_m": hist[-1][2], "iterations": hist}
