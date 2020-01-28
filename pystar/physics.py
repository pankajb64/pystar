from math import exp

from astropy.units.si import kg, m, s, W, K, N, Pa

from pystar.composition import CNO
from .constants import A_bf, A_ff, g_ff, A_es, A_Hm, a, m_H, A_pp, f_pp, A_CNO, \
    A_He, f_3a
from astropy.units import dimensionless_unscaled
from astropy.constants import k_B


def specific_heat_ratio(n_atoms=1):
    assert n_atoms == 1, "Only Monoatomic gas supported for now"

    # return the specific heat Cp/Cv for monoatomic gas
    # Carroll & Ostlie Eqn (10.80)
    return 5 / 3


def guillotine_over_gaunt_factor(rho: kg / m ** 3, X: dimensionless_unscaled,
                                 surface: bool) -> kg ** (1 / 5) / m ** (
        3 / 5):
    """
    Compute guillotine-over-gaunt (t/g_bf)
    Based on "Introduction to stellar atmospheres and interiors"
    by Eva Novotny, 1973 P. 469

    StatStar uses different leading multiplier in different scenarios,
    using 2.82 once and 0.708 another time. The above book uses 2.82
    which I've stuck to, since I couldn't figure out where 0.708 comes from.
    Maybe its the inverse of the 1.4 factor used in the same formula, but I'm
    not sure.

    Also set it to be a small value if we're close to the surface, as per
    the model.
    """

    if surface:
        return 0.01 * (kg ** (1 / 5) / m ** (3 / 5))
    else:
        return 2.82 * ((rho * (1 + X)) ** 0.2)


def bound_free_absorption_coefficient(rho: kg / m ** 3,
                                      X: dimensionless_unscaled,
                                      Z: dimensionless_unscaled,
                                      surface: bool = False):
    """
    Compute the bound-free absorption coefficient,
    the multiplier for Kramer's law form (rho * T^3.5)
    """
    tog_bf = guillotine_over_gaunt_factor(rho, X, surface)
    return (A_bf / tog_bf) * Z * (1 + X)


def bound_free_opacity(T: K, rho: kg / m ** 3, X: dimensionless_unscaled,
                       Z: dimensionless_unscaled,
                       surface: bool = False) -> m ** 2 / kg:
    """
    Compute bound-free-opacity kappa_bf
    Carroll & Ostlie Eqn (9.22)
    """
    coeff_bf = bound_free_absorption_coefficient(rho, X, Z, surface)
    return coeff_bf * rho / (T ** 3.5)


def free_free_absorption_coefficient(X: dimensionless_unscaled,
                                     Z: dimensionless_unscaled):
    """
    Compute the free-free absorption coefficient,
    the multiplier for Kramer's law form (rho * T^3.5)
    """
    return A_ff * g_ff * (1 - Z) * (1 + X)


def free_free_opacity(T: K, rho: kg / m ** 3, X: dimensionless_unscaled,
                      Z: dimensionless_unscaled) -> m ** 2 / kg:
    """
    Compute free-free opacity kappa_ff
    Carroll & Ostlie Eqn (9.23)
    """
    coeff_ff = free_free_absorption_coefficient(X, Z)
    return coeff_ff * rho / (T ** 3.5)


def electron_scattering_opacity(X: dimensionless_unscaled) -> m ** 2 / kg:
    """
    Compute electron scattering opacity kappa_es
    Carroll & Ostlie Eqn (9.27)
    """
    return A_es * (1 + X)


def h_minus_ion_opacity(T: K, rho: kg / m ** 3, X: dimensionless_unscaled,
                        Z: dimensionless_unscaled) -> m ** 2 / kg:
    """
    Compute H-minus ion opacity kappa_Hminus
    Carroll & Ostlie Eqn (9.28)
    The StatStar model uses a bound of (1e-10, 1e-5) for density, whereas
    the book mentions (1e-7, 1e-2). I'm going by what the book says.
    The book also mentions the check X ~ 0.7 but the model doesn't.
    I am adding a check of X between 0.67 and 0.73
    """
    # Zero by default unless the condition below is satisfied
    kappa_h_minus = 0

    if (3000 * K <= T <= 6000 * K) and (
            1e-7 * (kg / m ** 3) <= rho <= 1e-2 * (kg / m ** 3)) and (
            0.67 < X < 0.73) and (0.001 < Z < 0.03):
        kappa_h_minus = A_Hm * (Z / 0.02) * (rho ** 0.5) * (T ** 9)

    return kappa_h_minus


def opacity(T: K, rho: kg / m ** 3, X: dimensionless_unscaled,
            Z: dimensionless_unscaled, surface=False) -> m ** 2 / kg:
    """
    Compute the total Rosseland mean opacity, kappa_bar
    """
    return bound_free_opacity(T, rho, X, Z, surface) + \
           free_free_opacity(T, rho, X, Z) + \
           electron_scattering_opacity(X) + \
           h_minus_ion_opacity(T, rho, X, Z)


def pt_log_gradient(P: Pa, T: K, P_prev: Pa, T_prev: K) -> Pa / K:
    """
    Compute the log pressure gradient w.r.t log temperature, i.e d(lnP)/d(lnT)
    This is equivalent to T/P (dP/dT). T and P are taken to be the average
    temperature and pressure, though if we take the ration, the 2 in the
    denominator cancels for both T and P, so we can simply represent the ratio
    of averages by the ratio of sums. dP is approximated as P_prev - P and
    dT as T_prev - T

    A value higher than 99.9 is truncated, similar to StatStar
    """

    dlnPdlnT = ((T + T_prev) / (P + P_prev)) * ((P_prev - P) / (T_prev - T))
    return 99.9 if dlnPdlnT > 99.9 else dlnPdlnT


def density(T: K, P: Pa, mu: dimensionless_unscaled) -> kg / m ** 3:
    """
    Compute the density using ideal gas equation of state
    Carroll & Ostlie Eqn(10.11)
    """
    P_rad = (1 / 3) * a * (T ** 4)
    P_gas = P - P_rad

    # I don't know why we're doing this
    # but also don't know whats the alternative
    if P_gas <= 0:
        P_gas = P

    # -1 to indicate error computing rho
    rho = -1
    if T > 0 and P > 0:
        rho = (P_gas * mu * m_H) / (k_B * T)

    return rho


def pp_chain_rate(T: K, rho: kg / m ** 3, X: dimensionless_unscaled) -> W / kg:
    """
    Compute Energy generation rate per unit mass for Proton-Proton Chains
    The coefficients psi_pp and C_pp are from
    Hansen and Kawaler, Eqns(6.65, 6.73, 6.74)
    Reference - https://www.google.com/books/edition/Stellar_Interiors/
    GI3qBwAAQBAJ?hl=en&gbpv=1&pg=PA302&printsec=frontcover&bsq=pp%20chain

    The pp-chain equation is Carroll & Ostlie Eqn(10.46)
    """
    # Temperature scale - 10^6 K
    T6 = T * 1e-6

    # psi_pp is a  correction factor that accounts for the simultaneous
    # occurrence of PP I, PP II, and PP III chains
    # Factor inside the exponential to deal with astropy units
    e_psi = -49.98 * K**(1/3)
    psi_pp = 1 + 1.412e8 * (1 / X - 1) * exp(e_psi * (T6 ** (-1 / 3)))

    # C_pp involves higher order correction terms
    # Factors defined separately to handle astropy units
    e1 = 0.0123 * K ** (-1 / 3)
    e2 = 0.0109 * K ** (-2 / 3)
    e3 = 0.000938 * K ** -1
    C_pp = 1 + e1 * (T6 ** (1 / 3)) + e2 * (
            T6 ** (2 / 3)) + e3 * T6

    # Factor inside the exponential to deal with astropy units
    e_pp = -33.80 * K ** (1 / 3)

    # pp-chain energy generation rate
    epsilon_pp = A_pp * rho * (X ** 2) * f_pp * psi_pp * C_pp * (
            T6 ** (-2 / 3)) * exp(e_pp * (T6 ** (-1 / 3)))

    return epsilon_pp


def cno_cycle_rate(T: K, rho: kg / m ** 3, X: dimensionless_unscaled,
                   Z: dimensionless_unscaled) -> W / kg:
    """
    Compute energy generation rate per unit mass for Carbon Nitrogen Oxygen
    Cycle.
    Computation of factors is from Kippenhahn and Weigert Eqn(18.65)
    Reference - https://archive.org/details/
    StellarStructureAndEvolutionKippenhahnWeigert/page/n87/mode/2up
    Computation of epsilon_cno is from Carroll & Ostlie Eqn(10.58)
    """
    # Temperature scale - 10^6 K
    T6 = T * 1e-6

    # Total Mass Fraction of C, N and O
    XCNO = CNO(Z)

    # Higher order correction term
    # Factors defined separately to handle astropy units
    e1 = 0.0027 * K**(-1/3)
    e2 = 0.00778 * K**(-2/3)
    e3 = 0.000149 * K**-1
    CCNO = 1 + e1 * (T6 ** (1 / 3)) - e2 * (T6 ** (2 / 3)) - e3 * T6

    # Factor inside the exponential to deal with astropy units
    e_CNO = -152.28 * K**(1/3)

    # CNO cycle energy generation rate
    epsilon_CNO = A_CNO * rho * X * XCNO * CCNO * (T6 ** (-2 / 3)) * exp(
        e_CNO * (T6 ** (-1 / 3)))

    return epsilon_CNO


def he_burning_rate(T: K, rho: kg/m**3, Y: dimensionless_unscaled) -> W/kg:
    """
    Compute energy generation rate per unit for Helium burning phase
    Computation of factors is from Kippenhahn and Weigert Eqn(18.67)
    Reference - https://archive.org/details/
    StellarStructureAndEvolutionKippenhahnWeigert/page/n87/mode/2up
    Computation of epsilon_cno is from Carroll & Ostlie Eqn(10.62)
    """
    # Temperature scale - 10^8 K
    T8 = T * 1e-8

    # Factor inside the exponential to deal with astropy units
    e_He = -44.027*K

    epsilon_He = A_He * (rho**2) * (Y**3) * (T8**-3) * f_3a * exp(e_He/T8)

    return epsilon_He


def energy_generation_rate(T: K, rho: kg/m**3, X: dimensionless_unscaled, Y: dimensionless_unscaled, Z: dimensionless_unscaled) -> W/kg:
    """
    Compute the nuclear energy generation rate per unit mass.
    Combines the energy generation rate for pp-chain, CNO cycle and He-burning
    """
    return pp_chain_rate(T, rho, X) + \
           cno_cycle_rate(T, rho, X, Z) + \
           he_burning_rate(T, rho, Y)
