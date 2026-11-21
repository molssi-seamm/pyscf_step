#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Driver for one spin-orbit task, run with the PySCF environment's Python.

    python so_driver.py input.json -o results.json

The task's directory holds this script, aoc.py and so_correction.py (flat). The
input names the task:

* ``"atom"``: an atom with an open p shell (p1 or p5), by aoc-SCF, scalar and two
  component (see aoc.p_shell_atom);
* ``"molecule"``: a closed-shell molecule, scalar RHF and two-component GHF with
  the SO-ECP (see so_correction.molecule_pair).

and gives the structure, the basis and the ECP (with its spin-orbit column) in
PySCF's internal format, so nothing here reads a basis-set library. Errors are
written to the results as ``{"error": ...}`` and the exit code is 1.
"""

import argparse
import json
import os
import sys
import time
import traceback


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input")
    parser.add_argument("-o", "--output", default="results.json")
    args = parser.parse_args(argv)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    t0 = time.time()
    out = {}
    try:
        data = json.load(open(args.input))
        out = run(data)
        rc = 0
    except Exception as e:  # report, do not lose the message
        out["error"] = f"{type(e).__name__}: {e}"
        out["traceback"] = traceback.format_exc()
        rc = 1
    out["wall_seconds"] = time.time() - t0
    with open(args.output, "w") as fd:
        json.dump(out, fd, indent=1)
    return rc


def _ecp(ecp):
    """JSON turns the element's [nelec, shells] into lists already; keep it."""
    return {element: [int(v[0]), v[1]] for element, v in ecp.items()}


def run(data):
    from pyscf import gto, lib

    import pyscf

    threads = int(data.get("threads") or os.environ.get("OMP_NUM_THREADS") or 1)
    lib.num_threads(threads)
    mol = gto.M(
        atom=[[a[0], tuple(a[1])] for a in data["atoms"]],
        basis=data["basis"],
        ecp=_ecp(data.get("ecp", {})),
        charge=int(data.get("charge", 0)),
        spin=int(data.get("spin", 0)),
        unit="Angstrom",
        max_memory=int(data.get("memory_mb", 4000)),
        verbose=int(data.get("verbose", 3)),
        output=data.get("log", "pyscf.log"),
    )
    out = {
        "task": data["task"],
        "pyscf": pyscf.__version__,
        "threads": lib.num_threads(),
        "nao": int(mol.nao),
        "nelectron": int(mol.nelectron),
    }
    if data["task"] == "atom":
        import aoc

        out.update(aoc.p_shell_atom(mol, int(data["n_p"])))
    elif data["task"] == "molecule":
        import so_correction

        out.update(so_correction.molecule_pair(mol))
    else:
        raise ValueError(f"unknown task '{data['task']}'")
    return out


if __name__ == "__main__":
    sys.exit(main())
