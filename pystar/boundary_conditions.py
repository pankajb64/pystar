from astropy.units.si import kg, W, m, Pa, K
from astropy.units import dimensionless_unscaled
from .composition import mean_molecular_weight
from .physics import specific_heat_ratio, bound_free_absorption_coefficient, \
    free_free_absorption_coefficient, pt_log_gradient
from astropy.constants import G, k_B, c
from .constants import m_H, a
from math import pi
from .stellar_equations import radiation_criterion

def surface_zone_temperature_radiative(M_s: kg, R_s: m,
                                       mu: dimensionless_unscaled, r: m) -> K:
    """
    Compute Temperature at the surface zone, assuming the surface is radiative
    Carroll & Ostlie Eqn (L.2)
    """
    return G * M_s * (mu * m_H / (4.25 * k_B)) * (1 / r - 1 / R_s)


def surface_zone_pressure_radiative(M_s: kg, L_s: W,
                                    mu: dimensionless_unscaled,
                                    A, T: K) -> Pa:
    """
    Compute Pressure at the surface zone, assuming the surface is radiative
    Carroll & Ostlie Eqn (L.1)
    """
    p0 = (1 / 4.25) * (16 * pi / 3) * (G * M_s / L_s) * (
            a * c * k_B / (A * mu * m_H))

    return (p0 ** 0.5) * (T ** 4.25)


def surface_zone_temperature_convective(M_s: kg, R_s: m,
                                        mu: dimensionless_unscaled, r: m,
                                        gamma: dimensionless_unscaled) -> K:
    """
    Compute Temperature at the surface zone, assuming the surface is convective
    Carroll & Ostlie Eqn (L.3)
    """
    gamma_ratio = gamma / (gamma - 1)
    return (G * M_s / gamma_ratio) * (mu * m_H / k_B) * (1 / r - 1 / R_s)


def surface_zone_pressure_convective(T: K, T_prev: K, P_prev: Pa,
                                     gamma: dimensionless_unscaled) -> Pa:
    """
    Compute Pressure at the surface zone, assuming the surface is convective
    Constant of proportionality is computed from previous pressure
    and temperature
    Carroll & Ostlie Eqn (10.83)
    """
    gamma_ratio = gamma / (gamma - 1)
    k_p_adiabatic = P_prev / (T_prev ** gamma_ratio)
    return k_p_adiabatic * (T, gamma_ratio)


def surface_zone_conditions(M_s: kg, L_s: W, M_prev: kg, L_prev: W, r_prev: m,
                            P_prev: Pa, T_prev: K, X: dimensionless_unscaled,
                            Y: dimensionless_unscaled,
                            Z: dimensionless_unscaled, dr: m,
                            rho: kg / m ** 3):
    """
    Try to compute the physical quantities in the surface zone,
    assuming pressure, temperature and density at the surface are zero

    Compute iteratively until the changes in Mass and Luminosity are
    within a certain acceptable range, or until a specific number of
    iterations have passed.
    """

    # maximum fraction change in M_r and L_r
    max_fraction_change: float = 1e-8

    # max number of iterations to run for
    max_iterations = 50

    r = r_prev + dr
    mu = mean_molecular_weight(X, Y, Z)
    gamma = specific_heat_ratio()
    gamma_ratio = gamma / (gamma - 1)

    # To begin with, assume the surface is radiative
    radiative = True

    T = surface_zone_temperature_radiative(M_s, r_prev, mu, r)
    A = bound_free_absorption_coefficient(rho, X, Z, surface=True) + \
        free_free_absorption_coefficient(X, Z)
    P = surface_zone_pressure_radiative(M_s, L_s, mu, A, T)

    # TODO try TF auto-diff
    dlnPdlnT = pt_log_gradient(P, T, P_prev, T_prev)

    # If Radiation Criterion not satisfied, then
    # recompute P and T based on Convection Eqns
    if not radiation_criterion(dlnPdlnT, gamma):
        radiative = False
        T = surface_zone_temperature_convective(M_s, r_prev, mu, r, gamma)
        P = surface_zone_pressure_convective(T, T_prev, P_prev, gamma)

    rho =


