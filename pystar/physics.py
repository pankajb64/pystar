from astropy.units.si import kg, m, s, W, K, N, Pa
from .constants import A_bf, A_ff, g_ff, A_es, A_Hm
from astropy.units import dimensionless_unscaled


def specific_heat_ratio(n_atoms=1):
    assert n_atoms == 1, "Only Monoatomic gas supported for now"

    # return the specific heat Cp/Cv for monoatomic gas
    # Carroll & Ostlie Eqn (10.80)
    return 5 / 3


def guillotine_over_gaunt_factor(rho: kg / m ** 3, X: dimensionless_unscaled,
                                 surface: bool) -> kg**(1/5)/m**(3/5):
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
        return 0.01 * (kg**(1/5)/m**(3/5))
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

    if (3000*K <= T <= 6000*K) and (1e-7*(kg/m**3) <= rho <= 1e-2*(kg/m**3)) and (
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


def pt_log_gradient(P: Pa, T: K, P_prev: Pa, T_prev: K) -> Pa/K:
    """
    Compute the log pressure gradient w.r.t log temperature, i.e d(lnP)/d(lnT)
    This is equivalent to T/P (dP/dT). T and P are taken to be the average
    temperature and pressure, though if we take the ration, the 2 in the
    denominator cancels for both T and P, so we can simply represent the ratio
    of averages by the ratio of sums. dP is approximated as P_prev - P and
    dT as T_prev - T

    A value higher than 99.9 is truncated, similar to StatStar
    """

    dlnPdlnT = ((T + T_prev)/(P + P_prev)) * ((P_prev - P)/(T_prev - T))
    return 99.9 if dlnPdlnT > 99.9 else dlnPdlnT


# def density()