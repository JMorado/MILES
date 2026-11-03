"""Periodic table data used for building and drawing molecules.

Lengths are in angstrom and masses in atomic mass units.
"""

import itertools
from dataclasses import dataclass

import numpy as np

# Z, symbol, name, standard atomic weight, covalent radius (Cordero 2008),
# van der Waals radius, display colour.
_TABLE = """
1 H Hydrogen 1.008 0.31 1.20 FFFFFF
2 He Helium 4.0026 0.28 1.40 D9FFFF
3 Li Lithium 6.94 1.28 1.82 CC80FF
4 Be Beryllium 9.0122 0.96 1.53 C2FF00
5 B Boron 10.81 0.84 1.92 FFB5B5
6 C Carbon 12.011 0.76 1.70 909090
7 N Nitrogen 14.007 0.71 1.55 3050F8
8 O Oxygen 15.999 0.66 1.52 FF0D0D
9 F Fluorine 18.998 0.57 1.47 90E050
10 Ne Neon 20.180 0.58 1.54 B3E3F5
11 Na Sodium 22.990 1.66 2.27 AB5CF2
12 Mg Magnesium 24.305 1.41 1.73 8AFF00
13 Al Aluminium 26.982 1.21 1.84 BFA6A6
14 Si Silicon 28.085 1.11 2.10 F0C8A0
15 P Phosphorus 30.974 1.07 1.80 FF8000
16 S Sulfur 32.06 1.05 1.80 FFFF30
17 Cl Chlorine 35.45 1.02 1.75 1FF01F
18 Ar Argon 39.948 1.06 1.88 80D1E3
19 K Potassium 39.098 2.03 2.75 8F40D4
20 Ca Calcium 40.078 1.76 2.31 3DFF00
21 Sc Scandium 44.956 1.70 2.15 E6E6E6
22 Ti Titanium 47.867 1.60 2.11 BFC2C7
23 V Vanadium 50.942 1.53 2.07 A6A6AB
24 Cr Chromium 51.996 1.39 2.06 8A99C7
25 Mn Manganese 54.938 1.39 2.05 9C7AC7
26 Fe Iron 55.845 1.32 2.04 E06633
27 Co Cobalt 58.933 1.26 2.00 F090A0
28 Ni Nickel 58.693 1.24 1.97 50D050
29 Cu Copper 63.546 1.32 1.96 C88033
30 Zn Zinc 65.38 1.22 2.01 7D80B0
31 Ga Gallium 69.723 1.22 1.87 C28F8F
32 Ge Germanium 72.630 1.20 2.11 668F8F
33 As Arsenic 74.922 1.19 1.85 BD80E3
34 Se Selenium 78.971 1.20 1.90 FFA100
35 Br Bromine 79.904 1.20 1.85 A62929
36 Kr Krypton 83.798 1.16 2.02 5CB8D1
37 Rb Rubidium 85.468 2.20 3.03 702EB0
38 Sr Strontium 87.62 1.95 2.49 00FF00
39 Y Yttrium 88.906 1.90 2.32 94FFFF
40 Zr Zirconium 91.224 1.75 2.23 94E0E0
41 Nb Niobium 92.906 1.64 2.18 73C2C9
42 Mo Molybdenum 95.95 1.54 2.17 54B5B5
43 Tc Technetium 98.0 1.47 2.16 3B9E9E
44 Ru Ruthenium 101.07 1.46 2.13 248F8F
45 Rh Rhodium 102.91 1.42 2.10 0A7D8C
46 Pd Palladium 106.42 1.39 2.10 006985
47 Ag Silver 107.87 1.45 2.11 C0C0C0
48 Cd Cadmium 112.41 1.44 2.18 FFD98F
49 In Indium 114.82 1.42 1.93 A67573
50 Sn Tin 118.71 1.39 2.17 668080
51 Sb Antimony 121.76 1.39 2.06 9E63B5
52 Te Tellurium 127.60 1.38 2.06 D47A00
53 I Iodine 126.90 1.39 1.98 940094
54 Xe Xenon 131.29 1.40 2.16 429EB0
55 Cs Caesium 132.91 2.44 3.43 57178F
56 Ba Barium 137.33 2.15 2.68 00C900
57 La Lanthanum 138.91 2.07 2.43 70D4FF
58 Ce Cerium 140.12 2.04 2.42 FFFFC7
59 Pr Praseodymium 140.91 2.03 2.40 D9FFC7
60 Nd Neodymium 144.24 2.01 2.39 C7FFC7
61 Pm Promethium 145.0 1.99 2.38 A3FFC7
62 Sm Samarium 150.36 1.98 2.36 8FFFC7
63 Eu Europium 151.96 1.98 2.35 61FFC7
64 Gd Gadolinium 157.25 1.96 2.34 45FFC7
65 Tb Terbium 158.93 1.94 2.33 30FFC7
66 Dy Dysprosium 162.50 1.92 2.31 1FFFC7
67 Ho Holmium 164.93 1.92 2.30 00FF9C
68 Er Erbium 167.26 1.89 2.29 00E675
69 Tm Thulium 168.93 1.90 2.27 00D452
70 Yb Ytterbium 173.05 1.87 2.26 00BF38
71 Lu Lutetium 174.97 1.87 2.24 00AB24
72 Hf Hafnium 178.49 1.75 2.23 4DC2FF
73 Ta Tantalum 180.95 1.70 2.22 4DA6FF
74 W Tungsten 183.84 1.62 2.18 2194D6
75 Re Rhenium 186.21 1.51 2.16 267DAB
76 Os Osmium 190.23 1.44 2.16 266696
77 Ir Iridium 192.22 1.41 2.13 175487
78 Pt Platinum 195.08 1.36 2.13 D0D0E0
79 Au Gold 196.97 1.36 2.14 FFD123
80 Hg Mercury 200.59 1.32 2.23 B8B8D0
81 Tl Thallium 204.38 1.45 1.96 A6544D
82 Pb Lead 207.2 1.46 2.02 575961
83 Bi Bismuth 208.98 1.48 2.07 9E4FB5
84 Po Polonium 209.0 1.40 1.97 AB5C00
85 At Astatine 210.0 1.50 2.02 754F45
86 Rn Radon 222.0 1.50 2.20 428296
87 Fr Francium 223.0 2.60 3.48 420066
88 Ra Radium 226.0 2.21 2.83 007D00
89 Ac Actinium 227.0 2.15 2.47 70ABFA
90 Th Thorium 232.04 2.06 2.45 00BAFF
91 Pa Protactinium 231.04 2.00 2.43 00A1FF
92 U Uranium 238.03 1.96 2.41 008FFF
93 Np Neptunium 237.0 1.90 2.39 0080FF
94 Pu Plutonium 244.0 1.87 2.43 006BFF
95 Am Americium 243.0 1.80 2.44 545CF2
96 Cm Curium 247.0 1.69 2.45 785CE3
97 Bk Berkelium 247.0 1.68 2.44 8A4FE3
98 Cf Californium 251.0 1.68 2.45 A136D4
99 Es Einsteinium 252.0 1.65 2.45 B31FD4
100 Fm Fermium 257.0 1.67 2.45 B31FBA
101 Md Mendelevium 258.0 1.73 2.46 B30DA6
102 No Nobelium 259.0 1.76 2.46 BD0D87
103 Lr Lawrencium 266.0 1.61 2.46 C70066
"""


@dataclass(frozen=True, slots=True)
class Element:
    """One chemical element.

    Attributes
    ----------
    number : int
        Atomic number (0 for the dummy atom).
    symbol, name : str
        Symbol and English name.
    mass : float
        Standard atomic weight.
    covalent_radius, vdw_radius : float
        Radii in angstrom.
    color : str
        Display colour as ``"#RRGGBB"``.
    row, column : int
        Position in the 18-column periodic table; lanthanides and actinides
        are in rows 8 and 9.
    """

    number: int
    symbol: str
    name: str
    mass: float
    covalent_radius: float
    vdw_radius: float
    color: str
    row: int
    column: int

    @property
    def is_metal(self) -> bool:
        """Whether the element is a metal (metalloids count as non-metals)."""
        return self.number not in _NONMETALS

    @property
    def is_f_block(self) -> bool:
        """Whether the element is a lanthanide or actinide."""
        return self.row >= 8


_NONMETALS = {1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86}

# First atomic number of each stretch of the table, with its row and the
# column of that first element. Lanthanides and actinides go in rows 8 and 9.
_LAYOUT = [(1, 1, 1), (2, 1, 18), (3, 2, 1), (5, 2, 13), (11, 3, 1), (13, 3, 13), (19, 4, 1),
           (37, 5, 1), (55, 6, 1), (57, 8, 3), (71, 6, 3), (87, 7, 1), (89, 9, 3), (103, 7, 3)]


def _table_position(z: int) -> tuple[int, int]:
    """Row and column of element ``z`` in the 18-column layout."""
    first, row, column = next(entry for entry in reversed(_LAYOUT) if entry[0] <= z)
    return row, column + z - first


def _load() -> dict[str, Element]:
    table = {}
    for line in _TABLE.strip().splitlines():
        z, symbol, name, mass, r_cov, r_vdw, color = line.split()
        row, column = _table_position(int(z))
        table[symbol] = Element(int(z), symbol, name, float(mass), float(r_cov), float(r_vdw),
                                "#" + color, row, column)
    return table


ELEMENTS: dict[str, Element] = _load()
"""Elements by symbol."""

BY_NUMBER: dict[int, Element] = {e.number: e for e in ELEMENTS.values()}
"""Elements by atomic number."""

DUMMY = Element(0, "X", "Dummy", 0.0, 0.5, 1.0, "#FF66CC", 0, 0)
"""Placeholder for unknown symbols."""

BY_Z: list[Element] = [BY_NUMBER.get(z, DUMMY) for z in range(max(BY_NUMBER) + 1)]
BY_Z[0] = DUMMY
"""Elements indexed by atomic number, the dummy atom at 0 (and in any gap)."""

SYMBOLS = np.array([e.symbol for e in BY_Z], dtype=object)
"""Element symbol by atomic number, for whole columns: ``SYMBOLS[mol.numbers]``."""
MASSES = np.array([e.mass or 1.0 for e in BY_Z])
"""Atomic mass by atomic number (1 for the dummy atom)."""
COVALENT_RADII = np.array([e.covalent_radius for e in BY_Z])
"""Covalent radius in angstrom by atomic number."""
VDW_RADII = np.array([e.vdw_radius for e in BY_Z])
"""Van der Waals radius in angstrom by atomic number."""


def element(symbol: str) -> Element:
    """Return the element with a given symbol, or :data:`DUMMY` if there is none."""
    return ELEMENTS.get(symbol, DUMMY)


def normalize_symbol(text: str) -> str | None:
    """Turn a label such as ``"c"``, ``"CL"``, ``"C12"``, ``"Fe(2+)"`` or ``"6"`` into a symbol.

    Parameters
    ----------
    text : str
        Atom label or atomic number.

    Returns
    -------
    str or None
        Element symbol, ``"X"`` for dummy atoms, or None if nothing matches.
    """
    text = text.strip()
    if text.isdigit():
        e = BY_NUMBER.get(int(text))
        return e.symbol if e else None
    letters = "".join(itertools.takewhile(str.isalpha, text))
    if not letters:
        return None
    if len(letters) >= 2 and letters[:2].capitalize() in ELEMENTS:
        return letters[:2].capitalize()
    one = letters[0].upper()
    return one if one in ELEMENTS or one == "X" else None


_STANDARD_VALENCE = {"H": 1, "B": 3, "C": 4, "N": 3, "O": 2, "F": 1, "Si": 4, "P": 3, "S": 2, "Cl": 1,
                     "Br": 1, "I": 1, "Se": 2, "Te": 2, "As": 3, "Ge": 4, "Sn": 4, "Pb": 4, "Sb": 3,
                     "Al": 3, "Ga": 3, "In": 3}
_VALENCE_BY_COLUMN = {1: 1, 2: 2, 13: 3, 14: 4, 15: 3, 16: 2, 17: 1, 18: 0}
_MAX_VALENCE = {"H": 1, "B": 4, "C": 4, "N": 4, "O": 2, "F": 1, "Si": 4, "P": 5, "S": 6, "Cl": 1,
                "Br": 1, "I": 3, "Se": 6, "As": 5}


def standard_valence(symbol: str) -> int:
    """Return the most common valence, used when adding hydrogens to loaded structures."""
    if symbol in _STANDARD_VALENCE:
        return _STANDARD_VALENCE[symbol]
    e = element(symbol)
    if e.column in _VALENCE_BY_COLUMN and not e.is_f_block:
        return _VALENCE_BY_COLUMN[e.column]
    return 6


def max_valence(symbol: str) -> int:
    """Return an upper bound on the bond-order sum, used to check perceived bonds."""
    return _MAX_VALENCE.get(symbol, 8)


# Reference bond lengths for common pairs, keyed by sorted symbols, then by
# bond order (1.5 means aromatic).
_BOND_LENGTHS: dict[tuple[str, str], dict[float, float]] = {
    ("C", "H"): {1: 1.09},
    ("H", "N"): {1: 1.01},
    ("H", "O"): {1: 0.96},
    ("H", "S"): {1: 1.34},
    ("H", "P"): {1: 1.42},
    ("H", "Si"): {1: 1.48},
    ("B", "H"): {1: 1.19},
    ("H", "H"): {1: 0.74},
    ("C", "C"): {1: 1.54, 1.5: 1.40, 2: 1.34, 3: 1.20},
    ("C", "N"): {1: 1.47, 1.5: 1.34, 2: 1.28, 3: 1.16},
    ("C", "O"): {1: 1.43, 1.5: 1.36, 2: 1.21, 3: 1.13},
    ("C", "S"): {1: 1.82, 1.5: 1.71, 2: 1.60},
    ("C", "F"): {1: 1.35},
    ("C", "Cl"): {1: 1.77},
    ("Br", "C"): {1: 1.94},
    ("C", "I"): {1: 2.14},
    ("C", "P"): {1: 1.84, 2: 1.67},
    ("C", "Si"): {1: 1.87},
    ("B", "C"): {1: 1.57},
    ("N", "N"): {1: 1.45, 1.5: 1.33, 2: 1.25, 3: 1.10},
    ("N", "O"): {1: 1.40, 1.5: 1.22, 2: 1.21},
    ("O", "O"): {1: 1.48, 2: 1.21},
    ("O", "S"): {1: 1.57, 2: 1.43},
    ("O", "P"): {1: 1.60, 2: 1.48},
    ("N", "S"): {1: 1.68, 2: 1.54},
    ("S", "S"): {1: 2.05},
    ("O", "Si"): {1: 1.64},
}


def _order_factor(order: float) -> float:
    """Shortening of a bond relative to a single bond."""
    return 1.0 - 0.115 * (order - 1.0)


def bond_length(a: str, b: str, order: float = 1.0) -> float:
    """Return the typical length of a bond.

    Parameters
    ----------
    a, b : str
        Element symbols.
    order : float, optional
        Bond order (1.5 for aromatic).

    Returns
    -------
    float
        Tabulated length when known, otherwise the sum of covalent radii
        shortened for multiple bonds.
    """
    table = _BOND_LENGTHS.get((a, b) if a <= b else (b, a))
    if table:
        return table[order] if order in table else table[1] * _order_factor(order)
    return (element(a).covalent_radius + element(b).covalent_radius) * _order_factor(order)


REFERENCE_ORDERS: dict[tuple[str, str], tuple[float | None, float | None, float | None, float | None]] = {
    ("C", "C"): (1.54, 1.40, 1.34, 1.20),
    ("C", "N"): (1.47, 1.34, 1.28, 1.16),
    ("C", "O"): (1.43, None, 1.21, 1.13),
    ("C", "S"): (1.82, 1.71, 1.60, None),
    ("N", "N"): (1.45, 1.33, 1.25, 1.10),
    ("N", "O"): (1.40, 1.22, None, None),
    ("O", "S"): (1.57, None, 1.43, None),
    ("O", "P"): (1.60, None, 1.48, None),
    ("O", "O"): (1.48, None, 1.21, None),
}
"""Lengths used to guess bond orders from geometry: (single, aromatic, double, triple)."""
