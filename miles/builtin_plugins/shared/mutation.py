"""Mutating amino acid residues: replacing a side chain by another one.

The new side chain comes from the bundled amino acid templates. The template
is superposed on the residue's backbone (N, CA, C), its side chain replaces
the old one, and the side-chain torsions chi1 and chi2 are scanned for the
least crowded conformation.
"""

import itertools
import math

import numpy as np

from miles import (
    AMINO_ACIDS,
    Atom,
    Molecule,
    add_hydrogens,
    amino_acid_template,
    clash_score,
    kabsch,
    near_mask,
    rotation_matrix,
    unit_vector,
)

BACKBONE = frozenset({"N", "CA", "C", "O", "OXT"})
"""Names of the backbone atoms of an amino acid, which a mutation keeps."""
"""Heavy backbone atoms, kept when a residue is mutated."""

CHI_STEPS = np.radians(np.arange(0.0, 360.0, 30.0))
ENVIRONMENT = 8.0
"""Atoms this close to CA are considered when choosing the side-chain conformation."""


def _names(mol: Molecule, atoms: list[int]) -> dict[str, int]:
    """Atom name to index within a residue (first occurrence of each name)."""
    names: dict[str, int] = {}
    for i in atoms:
        names.setdefault(mol.atoms[i].name, i)
    return names


def _backbone(mol: Molecule, residue: list[int]) -> set[int]:
    """Backbone heavy atoms of a residue and the hydrogens on them."""
    heavy = {i for i in residue if mol.atoms[i].name in BACKBONE}
    return heavy | {h for i in heavy for h in mol.hydrogens_on(i) if h in residue}


def _template_side_chain(template: Molecule, code: str) -> list[int]:
    """Template atoms beyond CA: the side chain with its hydrogens."""
    names = _names(template, list(range(len(template))))
    if code == "GLY":
        return []
    blocked = {names["CA"], names["N"]}
    return sorted(template.component(names["CB"], blocked=blocked))


def _chi_movers(template: Molecule, side_chain: list[int], code: str) -> list[tuple[int, int, list[int]]]:
    """Rotatable side-chain bonds as (axis start, axis end, atoms that turn), in template indices."""
    if code in ("GLY", "PRO"):
        return []
    names = _names(template, list(range(len(template))))
    ca, cb = names["CA"], names["CB"]
    chain = set(side_chain)
    torsions = [(ca, cb, [i for i in side_chain if i != cb])]
    gamma = [j for j in template.neighbors(cb) if j in chain and template.atoms[j].symbol != "H"]
    if len(gamma) == 1:
        cg = gamma[0]
        beyond = template.component(cg, blocked={cb})
        movers = [i for i in side_chain if i in beyond and i != cg]
        if any(template.atoms[i].symbol != "H" for i in movers):
            torsions.append((cb, cg, movers))
    return torsions


def _choose_rotamer(points: np.ndarray, torsions: list[tuple[int, int, list[int]]], environment: np.ndarray
                    ) -> np.ndarray:
    """Turn the side chain about its chi bonds to the least crowded conformation."""
    if not torsions or not len(environment):
        return points
    best, best_score = points, math.inf
    for angles in itertools.product(CHI_STEPS, repeat=len(torsions)):
        trial = points.copy()
        for (a, b, movers), phi in zip(torsions, angles, strict=True):
            axis = unit_vector(trial[b] - trial[a])
            trial[movers] = (trial[movers] - trial[b]) @ rotation_matrix(axis, phi).T + trial[b]
        score = clash_score(trial, environment)
        if score < best_score - 1e-12:
            best, best_score = trial, score
    return best


def mutate_residue(mol: Molecule, i: int, code: str) -> list[int]:
    """Replace the side chain of the amino acid containing atom ``i``.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data, changed in place.
    i : int
        Any atom of the residue.
    code : str
        Three-letter code of the new amino acid, e.g. ``"TRP"``.

    Returns
    -------
    list of int
        Indices of the new side-chain atoms (empty for glycine).

    Raises
    ------
    ValueError
        If the atom is not in an amino acid with N, CA and C atoms.
    KeyError
        If there is no template for ``code``.
    """
    code = code.upper()
    template = amino_acid_template(code).mol
    resname = mol.atoms[i].resname
    residue = sorted(mol.residue_atoms(i))
    names = _names(mol, residue)
    if resname not in AMINO_ACIDS or not {"N", "CA", "C"} <= names.keys():
        raise ValueError("The atom is not part of an amino acid with a complete backbone")
    with_hydrogens = any(mol.atoms[j].symbol == "H" for j in residue)
    reference = mol.atoms[names["CA"]]

    # superpose the template on the residue's N, CA and C
    tnames = _names(template, list(range(len(template))))
    frame = [tnames["N"], tnames["CA"], tnames["C"]]
    P, Q = template.coords[frame], mol.coords[[names["N"], names["CA"], names["C"]]]
    R = kabsch(P - P.mean(axis=0), Q - Q.mean(axis=0))
    placed = (template.coords - P.mean(axis=0)) @ R.T + Q.mean(axis=0)
    side = [t for t in _template_side_chain(template, code)
            if with_hydrogens or template.atoms[t].symbol != "H"]

    # atoms to remove: the old side chain, a surplus hydrogen on CA (glycine) and on N (proline)
    doomed = set(residue) - _backbone(mol, residue)
    ca_hydrogens = mol.hydrogens_on(names["CA"])
    if code != "GLY" and len(ca_hydrogens) > 1:
        cb_position = placed[tnames["CB"]]
        doomed.add(min(ca_hydrogens, key=lambda h: float(np.linalg.norm(mol.coords[h] - cb_position))))
    if code == "PRO":
        n_hydrogens = mol.hydrogens_on(names["N"])
        if n_hydrogens:
            cd_position = placed[tnames["CD"]]
            doomed.add(min(n_hydrogens, key=lambda h: float(np.linalg.norm(mol.coords[h] - cd_position))))
    old_proline = resname == "PRO"
    mapping = mol.remove_atoms(doomed)
    n_atom, ca_atom = mapping[names["N"]], mapping[names["CA"]]
    kept = [mapping[j] for j in residue if mapping[j] >= 0]

    # add the new side chain
    new: dict[int, int] = {}
    for t in side:
        record = template.atoms[t]
        atom = Atom(record.symbol, record.geom, record.nsites, record.valence, record.inferred,
                    name=record.name, resname=code, resseq=reference.resseq, chain=reference.chain,
                    hetero=False, ss=reference.ss)
        new[t] = mol.add_atom(atom, placed[t])
    for (a, b), order in template.bonds.items():
        if a in new and b in new:
            mol.add_bond(new[a], new[b], order)
    if code != "GLY":
        mol.add_bond(ca_atom, new[tnames["CB"]], 1.0)
    if code == "PRO":
        mol.add_bond(n_atom, new[tnames["CD"]], 1.0)

    # turn chi1 and chi2 to the least crowded conformation
    added = [new[t] for t in side]
    row = {t: k for k, t in enumerate(side)} | {tnames["CA"]: len(side)}   # CA is the last row
    torsions = [(row[a], row[b], [row[t] for t in movers if t in row])
                for a, b, movers in _chi_movers(template, side, code)]
    if torsions:
        around = near_mask(mol.coords, mol.coords[[ca_atom]], ENVIRONMENT)
        around[added] = False
        around[ca_atom] = False
        points = np.vstack([mol.coords[added], mol.coords[[ca_atom]]])
        mol.set_positions(added, _choose_rotamer(points, torsions, mol.coords[around])[:len(side)])

    for j in kept:
        mol.atoms[j].resname = code
    if code != "GLY":
        for h in mol.hydrogens_on(ca_atom):
            mol.atoms[h].name = "HA"
    if with_hydrogens and code == "GLY":
        add_hydrogens(mol, [ca_atom])
        for name, h in zip(("HA2", "HA3"), mol.hydrogens_on(ca_atom), strict=False):
            mol.atoms[h].name = name
    if with_hydrogens and old_proline and code != "PRO":
        add_hydrogens(mol, [n_atom])
    return added
