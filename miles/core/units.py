"""Physical constants and unit conversions used across MILES.

MILES works in angstrom, kcal/mol, elementary charges and debye. Files
often use atomic units or electronvolts; these factors convert them. The
values are CODATA 2018.
"""

BOHR = 0.529177210903
"""One bohr in angstrom."""

HARTREE = 627.509474
"""One hartree in kcal/mol."""

HARTREE_EV = 27.211386245988
"""One hartree in electronvolts."""

ELECTRONVOLT = 23.060548
"""One electronvolt in kcal/mol."""

DEBYE_PER_AU = 2.541746473
"""One atomic unit of dipole moment (e times bohr) in debye."""

DEBYE_PER_E_ANGSTROM = 4.80320471
"""One elementary charge times one angstrom, in debye."""

COULOMB = 332.0637
"""Coulomb's constant in kcal/mol times angstrom per squared elementary charge."""
