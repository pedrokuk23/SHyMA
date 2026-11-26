"""Earth model, frames and ground-station geometry (SHyMA).

WGS-84 geodetic altitude by default (spherical Earth of the author's earlier ASMAT satellite tool
as an option), J2 gravity,
Earth rotation for the atmosphere-relative velocity, and the
elevation / slant-range geometry of the ASMAT satellite tool geometry module
re-used for tracking and telemetry stations.
"""
import math

MU = 398600.4418e9        # m^3/s^2
RE = 6378.137e3           # m
J2 = 1.08263e-3
OMEGA = 7.2921159e-5      # rad/s
C_LIGHT = 2.99792458e8
K_BOLTZ = 1.380649e-23   # exact SI value (2019)

STATIONS = {
    # name: (lat deg, lon deg, alt m)
    "CLA": (-2.3167, -44.3667, 50.0),     # Alcantara Launch Center
    "CLBI": (-5.9236, -35.1617, 30.0),    # Barreira do Inferno Launch Center
    "Woomera": (-30.95, 136.50, 150.0),   # HyShot launch range (reference only)
    "Esrange": (67.8933, 21.1067, 328.0), # Esrange Space Center (VSB-30 validation)
}


# ---------------------------------------------------------------------------
# Earth model selection. Default ('wgs84'): geodetic altitude on the WGS-84
# ellipsoid, J2 gravity, rotating. 'sphere': spherical Earth of equatorial
# radius (the convention of the ASMAT satellite tool). 'round': sphere of radius R2 with
# inverse-square gravity (NASA NESC round-Earth cases).
# ---------------------------------------------------------------------------
MODEL = "wgs84"
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)
R2_AUTHALIC = 6371007.1809      # sphere of equal area to WGS-84 (NESC round Earth)
_DEFAULTS = {"RE": 6378.137e3, "J2": 1.08263e-3, "OMEGA": 7.2921159e-5}


def set_model(name="wgs84", j2=True, rotation=True, radius=None):
    """Select the Earth model (module-wide)."""
    global MODEL, RE, J2, OMEGA
    MODEL = name
    RE = radius if radius is not None else (R2_AUTHALIC if name == "round" else _DEFAULTS["RE"])
    J2 = (1.08262998905e-3 if name == "wgs84" else _DEFAULTS["J2"]) if j2 else 0.0
    OMEGA = (7.292115e-5 if name == "wgs84" else _DEFAULTS["OMEGA"]) if rotation else 0.0


def reset_model():
    set_model("wgs84")


set_model("wgs84")


def altitude(x, y, z):
    """Altitude above the reference surface (geodetic for WGS-84)."""
    if MODEL != "wgs84":
        return math.sqrt(x * x + y * y + z * z) - RE
    return _geodetic(x, y, z)[2]


def _geodetic(x, y, z):
    p = math.hypot(x, y)
    lon = math.atan2(y, x)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    for _ in range(6):
        N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N if abs(lat) < 1.5 else z / math.sin(lat) - N * (1.0 - WGS84_E2)
        lat = math.atan2(z, p * (1.0 - WGS84_E2 * N / (N + h)))
    N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    h = p / math.cos(lat) - N if abs(lat) < 1.5 else z / math.sin(lat) - N * (1.0 - WGS84_E2)
    return lat, lon, h


def gravity(x, y, z):
    """Point mass + J2 acceleration (m/s^2) in an Earth-centred frame."""
    r2 = x * x + y * y + z * z
    r = math.sqrt(r2)
    k = -MU / (r2 * r)
    f = 1.5 * J2 * (RE / r) ** 2
    zr2 = 5.0 * z * z / r2
    ax = k * x * (1.0 + f * (1.0 - zr2))
    ay = k * y * (1.0 + f * (1.0 - zr2))
    az = k * z * (1.0 + f * (3.0 - zr2))
    return ax, ay, az


def llh_to_ecef(lat_deg, lon_deg, h):
    la, lo = math.radians(lat_deg), math.radians(lon_deg)
    if MODEL == "wgs84":
        N = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(la) ** 2)
        return ((N + h) * math.cos(la) * math.cos(lo), (N + h) * math.cos(la) * math.sin(lo),
                (N * (1.0 - WGS84_E2) + h) * math.sin(la))
    r = RE + h
    return (r * math.cos(la) * math.cos(lo), r * math.cos(la) * math.sin(lo), r * math.sin(la))


def eci_to_ecef(x, y, z, t):
    th = OMEGA * t
    c, s = math.cos(th), math.sin(th)
    return (c * x + s * y, -s * x + c * y, z)


def ecef_to_llh(x, y, z):
    if MODEL == "wgs84":
        la, lo, h = _geodetic(x, y, z)
        return math.degrees(la), math.degrees(lo), h
    r = math.sqrt(x * x + y * y + z * z)
    return math.degrees(math.asin(z / r)), math.degrees(math.atan2(y, x)), r - RE


def enu_basis(lat_deg, lon_deg):
    la, lo = math.radians(lat_deg), math.radians(lon_deg)
    e = (-math.sin(lo), math.cos(lo), 0.0)
    n = (-math.sin(la) * math.cos(lo), -math.sin(la) * math.sin(lo), math.cos(la))
    u = (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))
    return e, n, u


def look_angles(station, x, y, z):
    """Elevation (deg), azimuth (deg) and slant range (m) from a station to an ECEF point."""
    lat, lon, h = station
    sx, sy, sz = llh_to_ecef(lat, lon, h)
    dx, dy, dz = x - sx, y - sy, z - sz
    rng = math.sqrt(dx * dx + dy * dy + dz * dz)
    e, n, u = enu_basis(lat, lon)
    de = dx * e[0] + dy * e[1] + dz * e[2]
    dn = dx * n[0] + dy * n[1] + dz * n[2]
    du = dx * u[0] + dy * u[1] + dz * u[2]
    el = math.degrees(math.asin(du / rng))
    az = math.degrees(math.atan2(de, dn)) % 360.0
    return el, az, rng


def great_circle_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    h = math.sin((p2 - p1) / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2.0 * RE * math.asin(min(1.0, math.sqrt(h))) / 1000.0   # haversine


def link_margin_db(slant_range_m, f_hz=2.25e9, eirp_dbw=3.0, g_over_t_dbk=10.0,
                   rate_bps=2.0e6, ebn0_req_db=11.0, losses_db=6.0):
    """Telemetry link margin (dB) with the classical SMAD link equation
    (the same form used in the ASMAT satellite tool communication analysis)."""
    lam = C_LIGHT / f_hz
    fspl = 20.0 * math.log10(4.0 * math.pi * slant_range_m / lam)
    ebn0 = (eirp_dbw - fspl - losses_db + g_over_t_dbk
            - 10.0 * math.log10(K_BOLTZ) - 10.0 * math.log10(rate_bps))
    return ebn0 - ebn0_req_db
