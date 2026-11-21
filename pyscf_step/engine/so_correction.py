# -*- coding: utf-8 -*-

"""The molecular spin-orbit correction: two-component minus scalar SCF.

For a closed-shell molecule at a fixed geometry, with the same basis and ECP:

    E_SO(mol) = E(2c GHF, SO-ECP) - E(scalar RHF)

The core electrons left explicit by a small-core ECP (I: 4s4p4d) contribute a large
second-order term (~100 kJ/mol for I) that cancels against the same term in the
atoms. The molecular (second-order) correction to the atomization energy, with the
atomic J-averages taken from experiment as the composite schemes do, is therefore

    delta_SO = -E_SO(mol) + sum_atoms E_SO(J-avg, atom)

with the atoms by aoc-SCF (see :mod:`.aoc`). Bross et al. (J. Phys. Chem. A 127, 704
(2023)) keep this SCF level; it fits ATcT best.

Known offset (SEAMM 5877, report it with each result; do not tune it away): with
uncontracted cc-pVTZ-PP and experimental r_e these values sit 7-10 % below Bross's
DIRAC numbers (HI 2.19 vs 2.4, CH3I 2.33 vs 2.6, I2 7.57 vs 8.0 kJ/mol), likely from
the geometries (theirs CCSD(T)-F12b) and DIRAC's aoc and 2c-ECP implementation. The
basis is not the cause: HI is 2.190 / 2.189 / 2.199 with unc-cc-pVTZ-PP /
unc-aug-cc-pVTZ-PP / unc-cc-pVQZ-PP.
"""

import numpy as np

try:
    from .aoc import HARTREE_TO_KJ, _converge, _ghf_dm
except ImportError:  # run as flat files in a task's directory
    from aoc import HARTREE_TO_KJ, _converge, _ghf_dm


def molecule_pair(mol, check_scalar_limit=True):
    """Scalar RHF and two-component GHF (SO-ECP) energies of a closed-shell molecule.

    Returns {"E_scalar", "E_2c", "E_SO_kJ"} and, if asked, the scalar-limit check
    GHF(no SOC) - RHF, which must be ~0.
    """
    from pyscf import scf

    if mol.spin != 0:
        raise ValueError("the molecular spin-orbit correction is for closed shells")
    if not mol.has_ecp_soc():
        raise ValueError("the molecule's ECPs carry no spin-orbit terms")
    rhf = scf.RHF(mol)
    e_sc = _converge(rhf, None)
    ghf = scf.GHF(mol)
    ghf.with_soc = True
    e_2c = _converge(ghf, _ghf_dm(rhf.make_rdm1()))
    out = {
        "E_scalar": float(e_sc),
        "E_2c": float(e_2c),
        "E_SO_kJ": float((e_2c - e_sc) * HARTREE_TO_KJ),
    }
    if check_scalar_limit:
        g0 = scf.GHF(mol)
        out["GHF_noSOC_minus_RHF"] = float(
            _converge(g0, _ghf_dm(rhf.make_rdm1())) - e_sc
        )
    return out


def delta_so(e_so_molecule_kj, atom_terms_kj):
    """delta_SO = -E_SO(mol) + sum of the atoms' E_SO(J-avg), in kJ/mol.

    Parameters
    ----------
    e_so_molecule_kj : float
        E(2c) - E(scalar) of the molecule.
    atom_terms_kj : iterable of float
        E_SO(J-avg) of each atom with a spin-orbit ECP (once per atom).
    """
    return float(-e_so_molecule_kj + np.sum(list(atom_terms_kj)))
