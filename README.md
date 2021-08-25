# pystar
Model simulating a simple Zero-Age Main Sequence (ZAMS) Homogenous star, based on `StatStar` as provided in Appendix L of Carroll and Ostlie 2nd. Ed.

See 
- `stellar_equations.py` for implementation of eqns of stellar structure and evolution, such as Hydrostatic Equilibrium, Mass Continuity, Energy generation and transport.
- `surface_zone.py` for computing the conditions at the surface (code for core extrapolation to be added)
- `physics.py` for calculation of quantities such as opacity, density, reaction rates.

The recommended astropy units for each quantity are specified as type hints.
