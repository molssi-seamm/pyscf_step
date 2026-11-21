#!/usr/bin/env python
"""The spin-orbit ECP rebuilt from the j-weighted scalar data (no PySCF needed)."""

import pytest

from pyscf_step.engine import ecp_so


def scalar_ecp_from_published(element):
    """The scalar ECP in PySCF's format, as the BSE gives it."""
    ref = ecp_so.published(element)
    shells = [[-1, [[] for _ in range(7)]]]
    for lval, prims in sorted(ref["scalar"].items()):
        by_r = [[] for _ in range(7)]
        by_r[2] = [list(p) for p in prims]
        shells.append([lval, by_r])
    return [ref["core_electrons"], shells]


def test_published_iodine():
    ref = ecp_so.published("I")
    assert ref["ecp"] == "ECP28MDF"
    assert sorted(ref["so"]) == [1, 2, 3]
    assert ecp_so.published("Xx") is None


def test_rebuilt_so_matches_published():
    """-2c/l and +2c/(l+1) reproduce Stuttgart's SO block (SEAMM 5876: 3.8e-8)."""
    ecp, table = ecp_so.add_so_from_jweighted(scalar_ecp_from_published("I"))
    assert ecp_so.worst_difference(table, ecp_so.published("I")["so"]) < 1e-7
    p = [shell for shell in ecp[1] if shell[0] == 1][0][1][2]
    assert p[0] == pytest.approx([15.720141, 67.416239, -134.832478])
    assert p[1] == pytest.approx([15.208222, 134.807696, 134.807696])


def test_ecp_with_so_checks():
    ecp, info = ecp_so.ecp_with_so("I", scalar_ecp_from_published("I"))
    assert info["scalar"] == 0.0 and info["so"] < 1e-7
    assert all(len(p) == 3 for shell in ecp[1] if shell[0] > 0 for p in shell[1][2])


def test_refusals():
    with pytest.raises(ecp_so.ECPError, match="no spin-orbit data"):
        ecp_so.ecp_with_so("Xx", scalar_ecp_from_published("I"))
    bad = scalar_ecp_from_published("I")
    bad[1][2][1][2][0][1] *= 1.01  # change one scalar p coefficient
    with pytest.raises(ecp_so.ECPError, match="differs from the published"):
        ecp_so.ecp_with_so("I", bad)
    odd = scalar_ecp_from_published("I")
    odd[1][2][1][2].pop()
    with pytest.raises(ecp_so.ECPError, match="odd number"):
        ecp_so.add_so_from_jweighted(odd)
    wrong_core = scalar_ecp_from_published("I")
    wrong_core[0] = 46
    with pytest.raises(ecp_so.ECPError, match="replaces 46 electrons"):
        ecp_so.ecp_with_so("I", wrong_core)
