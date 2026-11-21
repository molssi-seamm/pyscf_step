#!/usr/bin/env python
"""Basis sets and ECPs from the BSE in PySCF's format (no PySCF needed)."""

import pytest

pytest.importorskip("basis_set_exchange")

from pyscf_step import basis  # noqa: E402
from pyscf_step.engine.ecp_so import ECPError  # noqa: E402
from pyscf_step.spin_orbit import n_functions  # noqa: E402


def test_hi_sizes_match_phase0():
    """Uncontracted cc-pVTZ-PP (I) and cc-pVTZ (H): 97 + 16 functions, as in 5877."""
    b, e, info = basis.basis_and_ecp(["H", "I"], "cc-pVTZ-PP", "cc-pVTZ")
    assert n_functions(b["I"]) == 97
    assert n_functions(b["H"]) == 16
    assert info["I"]["ecp"] == "ECP28MDF"
    assert info["I"]["so_check"]["so"] < 1e-7
    assert info["H"] == {"basis": "unc-cc-pVTZ"}
    assert set(e) == {"I"}


def test_ecp_has_so_column_and_local_part():
    _, e, _ = basis.basis_and_ecp(["I"], "cc-pVTZ-PP", "cc-pVTZ")
    nelec, shells = e["I"]
    assert nelec == 28
    assert [s[0] for s in shells] == [-1, 0, 1, 2, 3]
    for lval, by_r in shells:
        for prim in by_r[2]:
            assert len(prim) == (3 if lval > 0 else 2)


def test_contracted_basis():
    b, _, info = basis.basis_and_ecp(["I"], "cc-pVTZ-PP", "cc-pVTZ", uncontract=False)
    assert info["I"]["basis"] == "cc-pVTZ-PP"
    assert n_functions(b["I"]) < 97


def test_element_without_so_data_is_refused():
    """Sn has an ECP in cc-pVTZ-PP but no SO data here yet: refused, not run scalar."""
    with pytest.raises(ECPError, match="no spin-orbit data"):
        basis.basis_and_ecp(["Sn"], "cc-pVTZ-PP", "cc-pVTZ")


def test_p_shells():
    assert basis.P_SHELL["I"] == 5 and basis.P_SHELL["Tl"] == 1
    assert "Te" not in basis.P_SHELL


def test_bromine():
    """Br's ECP10MDF with its SO part, checked against Stuttgart's block."""
    _, e, info = basis.basis_and_ecp(["Br"], "cc-pVTZ-PP", "cc-pVTZ")
    assert info["Br"]["ecp"] == "ECP10MDF" and e["Br"][0] == 10
    assert info["Br"]["so_check"]["so"] < 1e-6
