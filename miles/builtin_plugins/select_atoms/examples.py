"""Examples the Select dialog offers: common functional groups as SMARTS patterns."""

FUNCTIONAL_GROUPS: list[tuple[str, str]] = [
    ("Carbonyl C=O", "[CX3]=[OX1]"),
    ("Carboxylic acid or carboxylate", "[CX3](=[OX1])[OX2H1,OX1-]"),
    ("Amide", "[NX3][CX3]=[OX1]"),
    ("Hydroxyl (alcohol, phenol)", "[OX2H][#6;!$([#6]=[O,S,N])]"),
    ("Amine nitrogen", "[NX3,NX4+;!$(N[C,S]=[O,S,N]);!$(N#*);!$(N=*);!$(Na)]"),
    ("Aromatic ring atoms", "a"),
    ("Halogens", "[F,Cl,Br,I]"),
    ("Hydrogen-bond donors", "[N,O;!H0]"),
    ("Hydrogen-bond acceptors",
     "[$([O;!+]),$([n;H0;+0]),$([N;X1,X2;+0]),$([NX3;+0;!$(N[C,S]=[O,S,N]);!$(Na);!$(N-[#6]=[#6,#7])])]"),
    ("Sulfur", "[#16]"),
]
"""Common patterns offered in the Select dialog: (name, SMARTS)."""
