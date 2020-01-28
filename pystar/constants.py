from astropy.constants import sigma_sb, c
from astropy.units.si import kg, m, K, s
from astropy.units import dimensionless_unscaled

# Radiation constant
a = 4 * sigma_sb / c

# Mass of Hydrogen atom, as per Carroll & Ostlie (2nd Ed.)
m_H = 1.673532499 * kg

# Bound-free absorption coefficient
A_bf = 4.34e21 * ((K**3.5 * m**(22/5))/kg**(9/5))

# Free-free absorption coefficient
A_ff = 3.68e18 * ((K**3.5 * m**5)/kg**2)

# Gaunt factor, free-free absorption
g_ff = 1 * dimensionless_unscaled

# Electron-scattering absorption coefficient
A_es = 0.02 * (m**2/kg)

# Hminus ion absorption coefficient
A_Hm = 7.9e-34 * m**3.5/(K**9 * kg**1.5)

# pp-chain screening factor
f_pp = 1

# pp-chain reaction rate coefficient
A_pp = 0.241 * (K**(2/3) * m**5) / (kg * s**3)

# CNO cycle reaction rate coefficient
A_CNO = 8.67e20 * (K**(2/3) * m**5) / (kg * s**3)

# He-burning reaction rate coefficient
A_He = 50.9 * (K**3 * m**8)/(kg**2 * s**3)

# Screening factor for triple alpha process
f_3a = 1
