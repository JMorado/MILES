"""Code that several built-in plugins share.

This package is not a plugin and not part of the public API: it holds the
chemistry that more than one built-in feature needs but that the rest of the
application does not use, so it stays out of :mod:`miles.core`. Like the
plugins, it imports only :mod:`miles` and :mod:`miles.plugin`; plugins import
it, never each other. Code that only one plugin uses lives in that plugin.

:mod:`~miles.builtin_plugins.shared.hbonds`
    What a hydrogen bond is (Hydrogen Bonds, Binding Pocket).
:mod:`~miles.builtin_plugins.shared.mutation`
    Replacing a residue's side chain (Mutate Residue, Prepare Protein,
    and Build from Sequence through its peptide builder).
:mod:`~miles.builtin_plugins.shared.molecules`
    Whether a structure is a small molecule and has all its hydrogens
    (Properties, 2D drawings).
"""
