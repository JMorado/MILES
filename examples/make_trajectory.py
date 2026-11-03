"""Write alanine_dipeptide_md.xyz: a short molecular dynamics run of alanine dipeptide.

Run from the repository root with the development environment:

    python examples/make_trajectory.py

Forces come from Open Babel's MMFF94; the integrator is velocity Verlet with
a Berendsen thermostat at 300 K and 0.5 fs steps. The file has 200 frames,
one every 20 fs (4 ps in total). It is meant for trying the timeline and the
movie export, not for analysis.
"""

from pathlib import Path

import numpy as np
from openbabel import openbabel as ob

import miles

SMILES = "CC(=O)N[C@@H](C)C(=O)NC"   # Ace-Ala-NMe
TEMPERATURE = 300.0   # K
TIMESTEP = 0.5        # fs
COUPLING = 100.0      # fs, Berendsen time constant
FRAMES = 200
STEPS_PER_FRAME = 40
BOLTZMANN = 0.0019872041   # kcal/mol/K
ACCELERATION = 4.184e-4    # (kcal/mol/A) / amu in A/fs^2
SEED = 7


def main() -> None:
    """Run the dynamics and write the trajectory next to this script."""
    mol, _ = miles.from_smiles(SMILES)
    miles.clean_geometry(mol)
    obmol = miles.to_openbabel(mol)
    ff = ob.OBForceField.FindForceField("MMFF94")
    if ff is None or not ff.Setup(obmol):
        raise SystemExit("MMFF94 is not available")
    atoms = list(ob.OBMolAtomIter(obmol))
    masses = mol.masses()[:, None]

    def accelerations(x: np.ndarray) -> np.ndarray:
        for atom, position in zip(atoms, x, strict=True):
            atom.SetVector(*position)
        ff.SetCoordinates(obmol)
        ff.Energy(True)
        forces = np.array([[g.GetX(), g.GetY(), g.GetZ()] for g in (ff.GetGradient(a) for a in atoms)])
        return forces / masses * ACCELERATION

    def temperature(v: np.ndarray) -> float:
        return float((masses * v * v).sum() / ACCELERATION / (3 * len(v) * BOLTZMANN))

    rng = np.random.default_rng(SEED)
    x = mol.coords.copy()
    v = rng.normal(size=x.shape) * np.sqrt(BOLTZMANN * TEMPERATURE / masses * ACCELERATION)
    v -= (masses * v).sum(axis=0) / masses.sum()
    a = accelerations(x)
    lines: list[str] = []
    for frame in range(FRAMES):
        for _ in range(STEPS_PER_FRAME):
            x = x + v * TIMESTEP + 0.5 * a * TIMESTEP**2
            a_new = accelerations(x)
            v = v + 0.5 * (a + a_new) * TIMESTEP
            a = a_new
            v *= np.sqrt(1 + TIMESTEP / COUPLING * (TEMPERATURE / max(temperature(v), 1.0) - 1))
        t = (frame + 1) * STEPS_PER_FRAME * TIMESTEP / 1000
        lines.append(str(len(mol)))
        lines.append(f"alanine dipeptide, MMFF94 at {TEMPERATURE:.0f} K, t = {t:.2f} ps")
        lines.extend(f"{atom.symbol:2s} {px:12.6f} {py:12.6f} {pz:12.6f}"
                     for atom, (px, py, pz) in zip(mol.atoms, x, strict=True))
    out = Path(__file__).with_name("alanine_dipeptide_md.xyz")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out}: {FRAMES} frames, final temperature {temperature(v):.0f} K")


if __name__ == "__main__":
    main()
