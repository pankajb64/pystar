from astropy.units import dimensionless_unscaled


def mean_molecular_weight(X: dimensionless_unscaled, Y: dimensionless_unscaled,
                          Z: dimensionless_unscaled) -> dimensionless_unscaled:
    """
    Calculate mean molecular weight of the gas, assuming complete ionization.
    Carroll & Ostlie Eqn (10.16)
    """
    return 1/(2*X + 3*Y/4 + Z/2)


def CNO(Z: dimensionless_unscaled) -> dimensionless_unscaled:
    """
    Compute the total mass fraction of Carbon, Nitrogen and Oxygen.
    Here we just approximate it as the mass fraction of all heavy elements
    """
    return Z/2
