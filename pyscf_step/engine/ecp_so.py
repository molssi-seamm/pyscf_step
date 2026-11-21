# -*- coding: utf-8 -*-

"""Spin-orbit parts of energy-consistent (Stuttgart-Cologne) pseudopotentials.

The basis-set libraries (the Basis Set Exchange, PySCF's own) store these ECPs as
scalar potentials only. For the MDF potentials, each l > 0 term is kept as two
primitives, j = l - 1/2 then j = l + 1/2, weighted by l/(2l+1) and (l+1)/(2l+1) so
that their sum is the averaged potential. The j-dependent potentials, and from them
the spin-orbit operator, can therefore be rebuilt from the scalar data.

PySCF takes the spin-orbit part as a third coefficient per primitive,
``[exponent, c_scalar, c_so]``, for the operator U_SO(r) l.s with

    U_SO = 2/(2l+1) [V(l, l+1/2) - V(l, l-1/2)],

which is -2c/l on the j = l - 1/2 primitive and +2c/(l+1) on the j = l + 1/2 one
(c being the stored, weighted scalar coefficient). This is the Molpro convention
for the published SO blocks, so the rebuilt values can be checked against them.

This module needs only the standard library.
"""

import json
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "so_ecp" / "stuttgart_so.json"


class ECPError(Exception):
    pass


def published(element):
    """The published SO and scalar blocks for an element, or None.

    Returns {"ecp", "core_electrons", "reference", "source", "scalar", "so"}, with
    the blocks as {l: [(exponent, coefficient), ...]}.
    """
    data = json.loads(DATA.read_text())
    entry = data.get(element)
    if entry is None:
        return None
    out = dict(entry)
    for key in ("scalar", "so"):
        out[key] = {
            int(lval): [tuple(p) for p in prims] for lval, prims in entry[key].items()
        }
    return out


def add_so_from_jweighted(ecp):
    """Add the spin-orbit column to an ECP in PySCF's format.

    Parameters
    ----------
    ecp : [int, [[l, [[[exp, c], ...] for each power of r]], ...]]
        A scalar ECP as PySCF holds it (``l = -1`` is the local part).

    Returns
    -------
    (ecp_with_so, so_table)
        The ECP with ``[exp, c, c_so]`` primitives for l > 0, and the SO terms as
        {l: [(exponent, coefficient), ...]} for checking.
    """
    nelec, shells = ecp
    new = []
    table = {}
    for shell in shells:
        lval, by_r = shell[0], shell[1]
        if lval < 1:
            new.append([lval, [[list(p) for p in prims] for prims in by_r]])
            continue
        new_by_r = []
        for power, prims in enumerate(by_r):
            if not prims:
                new_by_r.append([])
                continue
            if len(prims) % 2:
                raise ECPError(
                    f"The l={lval} r^{power} term has an odd number of primitives, "
                    "so it is not stored as j-weighted pairs"
                )
            out = []
            for k in range(0, len(prims), 2):
                (e_lo, c_lo), (e_hi, c_hi) = prims[k][:2], prims[k + 1][:2]
                so_lo = -2.0 * c_lo / lval
                so_hi = 2.0 * c_hi / (lval + 1)
                out.append([e_lo, c_lo, so_lo])
                out.append([e_hi, c_hi, so_hi])
                table.setdefault(lval, []).extend([(e_lo, so_lo), (e_hi, so_hi)])
            new_by_r.append(out)
        new.append([lval, new_by_r])
    return [nelec, new], table


def scalar_table(ecp, power=2):
    """The scalar terms of one power of r, as {l: [(exponent, coefficient), ...]}."""
    return {
        lval: [tuple(p[:2]) for p in by_r[power]]
        for lval, by_r in ecp[1]
        if lval >= 0 and len(by_r) > power and by_r[power]
    }


def worst_difference(ours, reference):
    """The largest difference between two {l: [(exponent, coefficient)]} tables:
    absolute in the exponents, relative (to max(1, |c|)) in the coefficients."""
    worst = 0.0
    for lval, prims in reference.items():
        mine = ours.get(lval, [])
        if len(mine) != len(prims):
            raise ECPError(f"l={lval} has {len(mine)} terms, not {len(prims)}")
        for (e1, c1), (e2, c2) in zip(sorted(mine), sorted(prims)):
            worst = max(worst, abs(e1 - e2), abs(c1 - c2) / max(1.0, abs(c2)))
    if set(ours) - set(reference):
        raise ECPError(f"extra angular momenta {sorted(set(ours) - set(reference))}")
    return worst


def ecp_with_so(element, scalar_ecp, tolerance=2e-6):
    """The ECP with its spin-orbit part, checked against the published blocks.

    Refuses an element with no published SO data, or whose scalar ECP is not the
    published one, rather than running without spin-orbit.

    Returns (ecp_with_so, {"scalar": worst, "so": worst, "reference": ...}).
    """
    ref = published(element)
    if ref is None:
        raise ECPError(
            f"There is no spin-orbit data for {element}'s pseudopotential, so a "
            "two-component calculation with it would silently miss spin-orbit"
        )
    if scalar_ecp[0] != ref["core_electrons"]:
        raise ECPError(
            f"{element}: the ECP replaces {scalar_ecp[0]} electrons, but the "
            f"spin-orbit data are for {ref['ecp']} ({ref['core_electrons']} electrons)"
        )
    worst_scalar = worst_difference(scalar_table(scalar_ecp), ref["scalar"])
    if worst_scalar > tolerance:
        raise ECPError(
            f"{element}: the scalar ECP differs from the published {ref['ecp']} "
            f"(worst difference {worst_scalar:.2e})"
        )
    ecp, table = add_so_from_jweighted(scalar_ecp)
    worst_so = worst_difference(table, ref["so"])
    if worst_so > tolerance:
        raise ECPError(
            f"{element}: the spin-orbit terms rebuilt from the ECP differ from the "
            f"published ones (worst difference {worst_so:.2e})"
        )
    return ecp, {"scalar": worst_scalar, "so": worst_so, "reference": ref["reference"]}
