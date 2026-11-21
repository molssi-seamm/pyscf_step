# -*- coding: utf-8 -*-

"""Basis sets and ECPs from the Basis Set Exchange, in PySCF's internal format.

This runs in the SEAMM environment (which has basis_set_exchange) and needs no
PySCF: the task's input then carries the basis and the ECP, with its spin-orbit
part, so the PySCF side reads no library. The spin-orbit terms are rebuilt and
checked by :mod:`pyscf_step.engine.ecp_so`.
"""

from .engine import ecp_so


class BasisError(Exception):
    pass


def _bse():
    import basis_set_exchange

    return basis_set_exchange


def has_element(name, symbol):
    """Whether the BSE basis ``name`` covers the element."""
    bse = _bse()
    try:
        bse.get_basis(name, elements=[symbol])
    except Exception:
        return False
    return True


def _element_data(name, symbol):
    bse = _bse()
    from basis_set_exchange import lut

    z = lut.element_Z_from_sym(symbol)
    data = bse.get_basis(name, elements=[symbol])
    return data["elements"][str(z)]


def electron_basis(name, symbol, uncontract=False):
    """The element's basis as PySCF's list [[l, [exp, c...], ...], ...]."""
    shells = []
    for shell in _element_data(name, symbol).get("electron_shells", []):
        ams = shell["angular_momentum"]
        exps = [float(x) for x in shell["exponents"]]
        coefs = [[float(x) for x in row] for row in shell["coefficients"]]
        if len(ams) == 1:
            columns = {ams[0]: coefs}
        else:  # an sp-type shell: one coefficient row per angular momentum
            columns = {lval: [coefs[i]] for i, lval in enumerate(ams)}
        for lval, rows in columns.items():
            shells.append(
                [lval] + [[e] + [row[i] for row in rows] for i, e in enumerate(exps)]
            )
    if uncontract:
        by_l = {}
        for shell in shells:
            for prim in shell[1:]:
                by_l.setdefault(shell[0], set()).add(prim[0])
        shells = [
            [lval, [e, 1.0]]
            for lval in sorted(by_l)
            for e in sorted(by_l[lval], reverse=True)
        ]
    return shells


def scalar_ecp(name, symbol):
    """The element's ECP as PySCF's [nelec, [[l, by_r], ...]], or None.

    The local part (stored by the BSE with the highest angular momentum, one above
    the others) becomes PySCF's l = -1. Only scalar potentials are accepted: the
    spin-orbit part is rebuilt separately and checked.
    """
    data = _element_data(name, symbol)
    pots = data.get("ecp_potentials")
    if not pots:
        return None
    if any(p.get("ecp_type", "scalar_ecp") != "scalar_ecp" for p in pots):
        raise BasisError(f"{name} for {symbol}: only scalar ECP terms are handled")
    lmax = max(p["angular_momentum"][0] for p in pots)
    others = sorted(
        p["angular_momentum"][0] for p in pots if p["angular_momentum"][0] != lmax
    )
    if others != list(range(lmax)):
        raise BasisError(f"{name} for {symbol}: unexpected ECP angular momenta")
    shells = []
    for p in pots:
        lval = p["angular_momentum"][0]
        by_r = [[] for _ in range(7)]
        for r, e, c in zip(
            p["r_exponents"], p["gaussian_exponents"], p["coefficients"][0]
        ):
            if float(c) != 0.0:
                by_r[int(r)].append([float(e), float(c)])
        shells.append([-1 if lval == lmax else lval, by_r])
    shells.sort(key=lambda s: s[0])
    return [int(data["ecp_electrons"]), shells]


def basis_and_ecp(symbols, basis, other_basis, uncontract=True, spin_orbit=True):
    """The basis and ECPs for a set of elements, as the driver takes them.

    Each element uses ``basis`` if the BSE has it there, else ``other_basis``. An
    element with an ECP gets the ECP's spin-orbit part, checked against the
    published data; an element with an ECP but no such data is refused when
    ``spin_orbit`` is true.

    Returns (basis {symbol: shells}, ecp {symbol: ecp}, info {symbol: {...}}).
    """
    out_basis, out_ecp, info = {}, {}, {}
    for symbol in sorted(set(symbols)):
        name = basis if has_element(basis, symbol) else other_basis
        if not has_element(name, symbol):
            raise BasisError(f"Neither {basis} nor {other_basis} covers {symbol}")
        out_basis[symbol] = electron_basis(name, symbol, uncontract=uncontract)
        entry = {"basis": ("unc-" if uncontract else "") + name}
        ecp = scalar_ecp(name, symbol)
        if ecp is not None:
            if spin_orbit:
                ecp, check = ecp_so.ecp_with_so(symbol, ecp)
                entry.update(
                    ecp=ecp_so.published(symbol)["ecp"],
                    so_check=check,
                    core_electrons=ecp[0],
                )
            out_ecp[symbol] = ecp
        info[symbol] = entry
    return out_basis, out_ecp, info


#: Valence p electrons of the neutral atoms whose aoc atom we can do (p1, p5).
P_SHELL = {
    "B": 1,
    "Al": 1,
    "Ga": 1,
    "In": 1,
    "Tl": 1,
    "F": 5,
    "Cl": 5,
    "Br": 5,
    "I": 5,
    "At": 5,
}
