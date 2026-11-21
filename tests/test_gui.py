# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialogs: the PySCF step and its Spin Orbit sub-step.
Skipped when no display is available."""

import pytest


def _root():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw

    Pmw.initialise(root)
    return root


def _tk_node(root, namespace, name):
    import seamm

    flowchart = seamm.Flowchart(namespace=namespace, directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace=namespace + ".tk"
    )
    node = flowchart.create_node(name)
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get(name)
    return plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )


def test_pyscf_dialog():
    root = _root()
    try:
        tk_node = _tk_node(root, "org.molssi.seamm", "PySCF")
        tk_node.create_dialog()
    finally:
        root.destroy()


def test_spin_orbit_dialog():
    root = _root()
    try:
        tk_node = _tk_node(root, "org.molssi.seamm.pyscf", "SpinOrbit")
        tk_node.create_dialog()
        tk_node.reset_dialog()
        for key in ("route", "method", "basis", "other basis", "uncontract"):
            assert tk_node[key].grid_info() != {}
        for value in ("cc-pVTZ-PP", "aug-cc-pVTZ-PP", "cc-pVQZ-PP"):
            tk_node["basis"].set(value)
            tk_node.reset_dialog()
    finally:
        root.destroy()
