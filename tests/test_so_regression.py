#!/usr/bin/env python
"""Regression against SEAMM 5877 (phase 0b, stamp 41d20942). Needs PySCF and runs
SCF calculations (about a minute on 8 cores), so it is skipped without PySCF and
marked slow: run it where quantum chemistry may run (CI or a cluster), not on a
laptop."""

import json
from pathlib import Path

import pytest

pyscf = pytest.importorskip("pyscf")

from pyscf import gto  # noqa: E402
from pyscf.gto.basis import parse_nwchem, parse_nwchem_ecp  # noqa: E402

from pyscf_step.engine import aoc, ecp_so, so_correction  # noqa: E402

DATA = Path(__file__).parent / "data" / "bse_ccpvtz_pp_I_HC.json"
pytestmark = pytest.mark.slow


def _mol(atoms, spin=0):
    bases = json.loads(DATA.read_text())
    ecp_text = "ECP" + bases["cc-pVTZ-PP"].split("ECP", 1)[1]
    ecp, _ = ecp_so.ecp_with_so("I", parse_nwchem_ecp.parse(ecp_text, "I"))
    basis = {
        "I": gto.uncontract(
            parse_nwchem.parse(bases["cc-pVTZ-PP"].split("ECP")[0], "I")
        ),
        "H": gto.uncontract(parse_nwchem.parse(bases["cc-pVTZ"], "H")),
    }
    return gto.M(
        atom=atoms,
        basis={k: basis[k] for k in {a[0] for a in atoms}},
        ecp={"I": ecp},
        spin=spin,
        unit="Angstrom",
        verbose=0,
    )


@pytest.fixture(scope="module")
def iodine():
    return aoc.p_shell_atom(_mol([["I", (0, 0, 0)]], spin=1), 5)


def test_iodine_atom(iodine):
    """5877: splitting 7560.43 cm-1, E_SO(J-avg) -103.048 kJ/mol; aoc checks pass."""
    assert iodine["splitting_cm"] == pytest.approx(7560.43, abs=0.5)
    assert iodine["E_SO_Javg_kJ"] == pytest.approx(-103.048, abs=0.005)
    for state in iodine["states"].values():
        assert state["determinant_average_check"] < 1e-8


def test_hi_delta_so(iodine):
    """5877: HI delta_SO 2.190 kJ/mol (Bross 2023 SCF 2.4; known offset -0.21)."""
    pair = so_correction.molecule_pair(_mol([["H", (0, 0, 0)], ["I", (0, 0, 1.6092)]]))
    assert abs(pair["GHF_noSOC_minus_RHF"]) < 1e-9
    assert -pair["E_SO_kJ"] == pytest.approx(105.239, abs=0.005)
    d = so_correction.delta_so(pair["E_SO_kJ"], [iodine["E_SO_Javg_kJ"]])
    assert d == pytest.approx(2.190, abs=0.005)
