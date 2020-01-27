from astropy.constants import G, c, k_B
from .constants import a, m_H
from math import pi
from astropy.units.si import kg, m, s, W, K, N, Pa
from astropy.units import dimensionless_unscaled


def hydrostatic_equilibrium(M_r: kg, rho: kg/m**3, r: m) -> Pa/m:
    """
    Hydrostatic Equilibrium Equation (Pressure Gradient)
    Carroll & Ostlie Eqn (10.6)
    """
    return - (G * M_r * rho) / r ** 2


def mass_continuity(r: m, rho: kg / m ** 3) -> kg / m:
    """
    Mass Continuity (or Mass Conservation) Equation (Mass Gradient)
    Carroll & Ostlie Eqn (10.7)
    """
    return 4 * pi * r ** 2 * rho


def energy_generation(r: m, rho: kg / m ** 3, epsilon: W / kg) -> W / m:
    """
    Energy generation Equation (Luminosity Gradient)
    Carroll & Ostlie Eqn (10.36)
    """
    return 4 * pi * r ** 2 * rho * epsilon


def radiation_transport(L_r: W, T: K, r: m, rho: kg / m ** 3,
                        kappa: m ** 2 / kg) -> K / m:
    """
    Radiation Transport Equation (Radiative Temperature Gradient)
    Carroll & Ostlie Eqn (10.68)
    """
    return (-3 * kappa * rho * L_r) / (4 * a * c * T ** 3 * 4 * pi * r ** 2)


def convection_transport(gamma: dimensionless_unscaled,
                         mu: dimensionless_unscaled, M_r: kg, r: m) -> K / m:
    """
    Convective or Adiabatic Transport Equation (Convective Temperature Gradient)
    Carroll & Ostlie Eqn (10.89)
    """
    gamma_ratio = gamma/(gamma - 1)
    return - (1/gamma_ratio) * (mu * m_H / k_B) * (G * M_r / r**2)


def radiation_criterion(dlnPdlnT: Pa/K,
                         gamma: dimensionless_unscaled) -> bool:
    """
    Criterion for radiation to occur
    Carroll & Ostlie Eqn (10.95)
    """
    gamma_ratio = gamma / (gamma - 1)
    return dlnPdlnT > gamma_ratio
