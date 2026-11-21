#!/usr/bin/env python
"""The Spin Orbit sub-step's parameters and description."""

from pyscf_step import SpinOrbit
from pyscf_step.spin_orbit import n_functions


def test_defaults():
    so = SpinOrbit()
    values = so.parameters.values_to_dict()
    assert values["route"] == "SO-ECP"
    assert values["method"] == "HF"
    assert values["basis"] == "cc-pVTZ-PP"
    assert values["other basis"] == "cc-pVTZ"
    assert values["uncontract"] == "yes"


def test_description():
    so = SpinOrbit()
    so._id = (1, 2)
    text = so.description_text()
    assert text.startswith("Step 1.2: Spin Orbit")
    assert "average-of-configuration" in text


def test_n_functions():
    assert (
        n_functions([[0, [1.0, 1.0]], [1, [2.0, 0.5, 0.5]], [2, [1.0, 1.0]]])
        == 1 + 6 + 5
    )
