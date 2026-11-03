"""Superposing one structure onto another and measuring the RMSD.

Atoms are paired first, then the Kabsch algorithm gives the rotation that
best overlays the pairs. Proteins are paired by aligning their sequences and
matching CA atoms; identical small molecules by their atom order or, when the
atoms are numbered differently, by their bond graphs.
"""

import difflib
from collections import Counter
from dataclasses import dataclass

import numpy as np

from miles import AMINO_ACIDS, ONE_LETTER, Molecule, kabsch

MIN_PAIRS = 3
REFINE_ROUNDS = 3
BRANCH_DEPTH = 3     # symmetric choices for which every partner is tried
MAX_MAPPINGS = 256   # candidate pairings compared by RMSD

Pairs = list[tuple[int, int]]


@dataclass(frozen=True)
class Superposition:
    """A rigid motion that overlays one structure on another.

    Attributes
    ----------
    rotation : numpy.ndarray
        Rotation matrix, shape ``(3, 3)``.
    translation : numpy.ndarray
        Shift applied after the rotation: ``x' = x @ rotation.T + translation``.
    rmsd : float
        Root-mean-square deviation of the paired atoms after the motion, in angstrom.
    pairs : list of tuple of int
        Paired atoms as ``(mobile index, target index)``.
    method : str
        How the atoms were paired.
    """

    rotation: np.ndarray
    translation: np.ndarray
    rmsd: float
    pairs: Pairs
    method: str


def _fit(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Rotation, translation and RMSD of the best overlay of points ``P`` on ``Q``."""
    p0, q0 = P.mean(axis=0), Q.mean(axis=0)
    R = kabsch(P - p0, Q - q0)
    t = q0 - p0 @ R.T
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(axis=1).mean()))
    return R, t, rmsd


# -- pairing proteins by sequence ----------------------------------------------------

def _alpha_carbons(mol: Molecule) -> tuple[list[int], str]:
    """CA atoms of amino acid residues in file order, with the one-letter sequence."""
    indices, letters = [], []
    seen: set[tuple[str, int, str]] = set()
    for i, a in enumerate(mol.atoms):
        key = (a.chain, a.resseq, a.resname)
        if a.name == "CA" and a.resname in AMINO_ACIDS and key not in seen:
            seen.add(key)
            indices.append(i)
            letters.append(ONE_LETTER.get(a.resname, "X"))
    return indices, "".join(letters)


def _pairs_by_sequence(mobile: Molecule, target: Molecule) -> Pairs:
    """CA atoms of residues matched by aligning the two sequences."""
    ca_m, seq_m = _alpha_carbons(mobile)
    ca_t, seq_t = _alpha_carbons(target)
    matcher = difflib.SequenceMatcher(None, seq_m, seq_t, autojunk=False)
    return [(ca_m[a + k], ca_t[b + k]) for a, b, size in matcher.get_matching_blocks() for k in range(size)]


def _pair_rmsd(mobile: Molecule, target: Molecule, pairs: Pairs) -> float:
    return _fit(mobile.coords[[m for m, _ in pairs]], target.coords[[t for _, t in pairs]])[2]


# -- pairing identical molecules by their graphs ---------------------------------------

def _heavy(mol: Molecule) -> list[int]:
    heavy = [i for i, a in enumerate(mol.atoms) if a.symbol != "H"]
    return heavy if len(heavy) >= MIN_PAIRS else list(range(len(mol)))


def _refine(labels: list[int], neighbors: list[list[int]]) -> list[int]:
    """Refine atom classes by their neighbours' classes until nothing changes (Morgan-like)."""
    count = len(set(labels))
    while True:
        signatures = [(labels[v], tuple(sorted(labels[w] for w in neighbors[v]))) for v in range(len(labels))]
        table = {s: k for k, s in enumerate(sorted(set(signatures)))}
        labels = [table[s] for s in signatures]
        if len(table) == count:
            return labels
        count = len(table)


def _individualize(labels: list[int], neighbors: list[list[int]], n: int, first: int,
                   second: int) -> list[int]:
    """Give two atoms (one of each molecule) a label of their own and refine again."""
    labels = labels[:]
    labels[first] = labels[second] = max(labels) + 1
    labels = _refine(labels, neighbors)
    if Counter(labels[:n]) != Counter(labels[n:]):
        raise ValueError("The structures are different molecules, so their atoms cannot be paired")
    return labels


def _ambiguous(labels: list[int]) -> list[int] | None:
    """Return the smallest class with more than one atom per molecule, if any."""
    members: dict[int, list[int]] = {}
    for v, label in enumerate(labels):
        members.setdefault(label, []).append(v)
    ambiguous = [vs for vs in members.values() if len(vs) > 2]
    return min(ambiguous, key=len) if ambiguous else None


def _isomorphisms(labels: list[int], neighbors: list[list[int]], n: int) -> list[list[int]]:
    """Complete labellings that pair every atom, branching over the first few symmetric choices.

    Equivalent atoms are told apart by individualising one pair at a time, as
    in canonical labelling. For the first :data:`BRANCH_DEPTH` choices every
    candidate partner is tried, so that the geometrically best of the
    symmetry-equivalent pairings can be picked; later choices take the first
    candidate.
    """
    results: list[list[int]] = []
    pending = [(labels, 0)]
    while pending and len(results) < MAX_MAPPINGS:
        labels, depth = pending.pop()
        group = _ambiguous(labels)
        if group is None:
            results.append(labels)
            continue
        first = next(v for v in group if v < n)
        candidates = [v for v in group if v >= n]
        if depth >= BRANCH_DEPTH:
            while group is not None:
                first = next(v for v in group if v < n)
                labels = _individualize(labels, neighbors, n, first, next(v for v in group if v >= n))
                group = _ambiguous(labels)
            results.append(labels)
            continue
        pending.extend((_individualize(labels, neighbors, n, first, t), depth + 1)
                       for t in reversed(candidates))
    return results


def _graph_pairs(mobile: Molecule, target: Molecule) -> tuple[Pairs, list[list[int]]]:
    """Pair the heavy atoms of two identical molecules by graph isomorphism.

    The two graphs are refined together, so equal labels mean equivalent
    atoms. Of the pairings found, the one with the lowest RMSD is kept.

    Returns
    -------
    pairs : list of tuple of int
        Paired atom indices.
    classes : list of list of int
        Pair indices grouped by their symmetry class; any reassignment within
        a class is also chemically valid.
    """
    atoms_m, atoms_t = _heavy(mobile), _heavy(target)
    if len(atoms_m) != len(atoms_t):
        raise ValueError("The structures have different numbers of atoms, so they cannot be matched")
    nodes = [(mobile, i) for i in atoms_m] + [(target, i) for i in atoms_t]
    index = {(id(mol), i): k for k, (mol, i) in enumerate(nodes)}
    neighbors = [[index[(id(mol), j)] for j in mol.neighbors(i) if (id(mol), j) in index] for mol, i in nodes]
    invariants = [(mol.atoms[i].symbol, len(mol.neighbors(i)), len(mol.hydrogens_on(i))) for mol, i in nodes]
    table = {s: k for k, s in enumerate(sorted(set(invariants)))}
    n = len(atoms_m)
    labels = _refine([table[s] for s in invariants], neighbors)
    if Counter(labels[:n]) != Counter(labels[n:]):
        raise ValueError("The structures are different molecules, so their atoms cannot be paired")

    def pairing(complete: list[int]) -> Pairs:
        partner = {complete[v]: v for v in range(n, 2 * n)}
        return [(atoms_m[v], atoms_t[partner[complete[v]] - n]) for v in range(n)]

    best = min((pairing(c) for c in _isomorphisms(labels, neighbors, n)),
               key=lambda pairs: _pair_rmsd(mobile, target, pairs))
    by_class: dict[int, list[int]] = {}
    for k in range(n):
        by_class.setdefault(labels[k], []).append(k)
    return best, list(by_class.values())


def _improve_by_geometry(mobile: Molecule, target: Molecule, pairs: Pairs, classes: list[list[int]]) -> Pairs:
    """Swap symmetry-equivalent partners to the nearest atoms after overlaying, if that lowers the RMSD."""
    best, best_rmsd = pairs, _pair_rmsd(mobile, target, pairs)
    current = list(pairs)
    for _ in range(REFINE_ROUNDS):
        R, t, _ = _fit(mobile.coords[[m for m, _ in current]], target.coords[[t for _, t in current]])
        moved = mobile.coords @ R.T + t
        for group in classes:
            if len(group) < 2:
                continue
            free = [current[k][1] for k in group]
            for k in sorted(group, key=lambda k: current[k][0]):
                m = current[k][0]
                nearest = min(free, key=lambda j: float(np.linalg.norm(moved[m] - target.coords[j])))
                free.remove(nearest)
                current[k] = (m, nearest)
        rmsd = _pair_rmsd(mobile, target, current)
        if rmsd < best_rmsd - 1e-9:
            best, best_rmsd = list(current), rmsd
    return best


# -- public API ------------------------------------------------------------------------

def match_atoms(mobile: Molecule, target: Molecule) -> tuple[Pairs, str]:
    """Pair the atoms of two structures for superposition.

    Parameters
    ----------
    mobile, target : Molecule
        The structure to move and the reference.

    Returns
    -------
    pairs : list of tuple of int
        ``(mobile index, target index)`` pairs.
    method : str
        How they were paired.

    Raises
    ------
    ValueError
        If fewer than three atoms can be paired.
    """
    if _alpha_carbons(mobile)[0] and _alpha_carbons(target)[0]:
        pairs = _pairs_by_sequence(mobile, target)
        if len(pairs) < MIN_PAIRS:
            raise ValueError("The protein sequences have fewer than three residues in common")
        return pairs, "CA atoms by sequence"
    if len(mobile) == len(target) and mobile.symbols == target.symbols:
        return [(i, i) for i in range(len(mobile))], "all atoms in order"
    pairs, classes = _graph_pairs(mobile, target)
    return _improve_by_geometry(mobile, target, pairs, classes), "heavy atoms by element and graph"


def superpose(mobile: Molecule, target: Molecule, pairs: Pairs | None = None) -> Superposition:
    """Find the rigid motion that best overlays ``mobile`` on ``target``.

    Parameters
    ----------
    mobile, target : Molecule
        The structure to move and the reference (neither is changed).
    pairs : list of tuple of int, optional
        Atoms to overlay; found with :func:`match_atoms` when omitted.

    Returns
    -------
    Superposition
        The motion and its RMSD.

    Raises
    ------
    ValueError
        If fewer than three atoms can be paired.
    """
    method = "given pairs"
    if pairs is None:
        pairs, method = match_atoms(mobile, target)
    if len(pairs) < MIN_PAIRS:
        raise ValueError("At least three pairs of atoms are needed to superpose")
    P = mobile.coords[[m for m, _ in pairs]]
    Q = target.coords[[t for _, t in pairs]]
    R, t, rmsd = _fit(P, Q)
    return Superposition(R, t, rmsd, list(pairs), method)


def apply(mol: Molecule, superposition: Superposition) -> None:
    """Move a structure by a superposition, in place."""
    mol.coords = mol.coords @ superposition.rotation.T + superposition.translation
