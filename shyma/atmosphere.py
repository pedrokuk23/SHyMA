"""Atmosphere model for SHyMA.

0-86 km : U.S. Standard Atmosphere 1976 (geopotential layers).
>86 km  : log-linear interpolation of the upper-atmosphere density table of the
          author's earlier satellite tool ASMAT (de Albuquerque, ITA M.Sc. dissertation, 2019)
          (minimum / mean / maximum solar activity), inherited from the
          original orbital module, with a smooth bridge between 86 and 100 km.
"""
import math

R_AIR = 287.05287      # J/(kg K)
GAMMA = 1.4
G0 = 9.80665
R0_GEOPOT = 6356.766   # km

# US76 base values (geopotential km, K/km, K, Pa)
_HB = [0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.8520]
_LB = [-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0]
_TB = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
_PB = [101325.0, 22632.06, 5474.889, 868.0187, 110.9063, 66.93887, 3.956420]

# Upper-atmosphere density table (from the ASMAT satellite tool, Constants.m): altitude (km), density (kg/m^3)
UPPER_ALT = [0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650,
             700, 750, 800, 850, 900, 950, 1000, 1250, 1500]
UPPER_RHO = {
    "min": [1.2, 4.61e-7, 1.65e-9, 1.78e-10, 3.35e-11, 8.19e-12, 2.34e-12, 7.32e-13,
            2.47e-13, 8.98e-14, 3.63e-14, 1.68e-14, 9.14e-15, 5.74e-15, 3.99e-15,
            2.96e-15, 2.28e-15, 1.8e-15, 1.44e-15, 1.17e-15, 4.67e-16, 2.3e-16],
    "mean": [1.2, 4.79e-7, 1.81e-9, 2.53e-10, 6.24e-11, 1.95e-11, 6.98e-12, 2.72e-12,
             1.13e-12, 4.89e-13, 2.21e-13, 1.04e-13, 5.15e-14, 2.72e-14, 1.55e-14,
             9.63e-15, 6.47e-15, 4.66e-15, 3.54e-15, 2.79e-15, 1.11e-15, 5.21e-16],
    "max": [1.2, 5.1e-7, 2.04e-9, 3.52e-10, 1.06e-10, 3.96e-11, 1.67e-11, 7.55e-12,
            3.61e-12, 1.8e-12, 9.25e-13, 4.89e-13, 2.64e-13, 1.47e-13, 8.37e-14,
            4.39e-14, 3.00e-14, 1.91e-14, 1.27e-14, 8.84e-15, 2.59e-15, 1.22e-15],
}
T_UPPER = 186.87  # K, used above 86 km only to define a nominal sound speed


def _us76(z_km):
    """Return T [K], p [Pa], rho [kg/m^3] for geometric altitude z <= 86 km."""
    h = R0_GEOPOT * z_km / (R0_GEOPOT + z_km)
    i = 0
    for k in range(len(_HB) - 1):
        if h >= _HB[k]:
            i = k
    if i > 6:
        i = 6
    dh = h - _HB[i]
    L = _LB[i]
    T = _TB[i] + L * dh
    if L == 0.0:
        p = _PB[i] * math.exp(-G0 * dh * 1000.0 / (R_AIR * _TB[i]))
    else:
        p = _PB[i] * (_TB[i] / T) ** (G0 * 1000.0 / (R_AIR * L))
    rho = p / (R_AIR * T)
    return T, p, rho


_T86, _P86, _RHO86 = _us76(86.0)


def _upper_rho(z_km, level):
    tab = UPPER_RHO[level]
    if z_km >= UPPER_ALT[-1]:
        return tab[-1]
    for k in range(1, len(UPPER_ALT)):
        if z_km < UPPER_ALT[k]:
            z0, z1 = UPPER_ALT[k - 1], UPPER_ALT[k]
            r0, r1 = tab[k - 1], tab[k]
            return math.exp(math.log(r0) + (math.log(r1) - math.log(r0)) * (z_km - z0) / (z1 - z0))
    return tab[-1]


def atmosphere(z_km, level="mean", rho_factor=1.0):
    """Return (T, p, rho, a) at geometric altitude z_km.

    level      : 'min' | 'mean' | 'max' (affects z > 86 km only)
    rho_factor : multiplicative density dispersion (Monte Carlo)
    """
    if z_km < 0.0:
        z_km = 0.0
    if z_km <= 86.0:
        T, p, rho = _us76(z_km)
    elif z_km < 100.0:
        r100 = _upper_rho(100.0, level)
        rho = math.exp(math.log(_RHO86) + (math.log(r100) - math.log(_RHO86)) * (z_km - 86.0) / 14.0)
        T = T_UPPER
        p = rho * R_AIR * T
    else:
        rho = _upper_rho(z_km, level)
        T = T_UPPER
        p = rho * R_AIR * T
    rho *= rho_factor
    p *= rho_factor
    a = math.sqrt(GAMMA * R_AIR * T)
    return T, p, rho, a


# ---------------------------------------------------------------------------
# Optional site- and date-specific atmosphere: NRLMSISE-00 (via pymsis).
# A table from 0 to 300 km is built once and interpolated (log-density).
# ---------------------------------------------------------------------------
def msis_atmosphere(lat, lon, date, f107=150.0, f107a=150.0, ap=4.0, dz=0.5):
    """Return an atmosphere callable (alt_km, level, rho_factor) based on
    NRLMSISE-00 at a fixed site and date (numpy-datetime string)."""
    import numpy as np
    from pymsis import msis
    alts = np.arange(0.0, 300.0 + dz, dz)
    out = np.squeeze(msis.run(np.datetime64(date), lon, lat, alts, f107s=[f107], f107as=[f107a],
                              aps=[[ap] * 7], version=0))
    rho = out[:, 0]
    temp = out[:, 10]
    ndens = np.nansum(out[:, 1:10], axis=1)          # number densities (1/m^3)
    pres = ndens * 1.380649e-23 * temp
    lrho = [math.log(float(v)) for v in rho]
    tab_T = [float(v) for v in temp]
    tab_p = [float(v) for v in pres]
    tab_rho = [float(v) for v in rho]

    def f(alt_km, level="mean", rho_factor=1.0):
        z = min(max(alt_km, 0.0), 300.0)
        i = min(int(z / dz), len(lrho) - 2)
        w = (z - i * dz) / dz
        r = math.exp(lrho[i] + w * (lrho[i + 1] - lrho[i])) * rho_factor
        T = tab_T[i] + w * (tab_T[i + 1] - tab_T[i])
        p = (tab_p[i] + w * (tab_p[i + 1] - tab_p[i])) * rho_factor
        a = math.sqrt(GAMMA * p / r) if r > 0 else math.sqrt(GAMMA * R_AIR * T)
        return T, p, r, a
    f.table = {"alt_km": [float(v) for v in alts], "rho": tab_rho, "T": tab_T, "p": tab_p}
    return f
