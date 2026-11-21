#!/usr/bin/env python
"""The Spin Orbit sub-step end to end: a flowchart run as a job runs it, through
seamm_exec, with PySCF in this Python (``pyscf.ini``: installation = local).

Slow (about a minute) and needs PySCF: it runs in CI (``.[test]``), not on a laptop.
"""

import re
import sys
import textwrap
from pathlib import Path

import pytest

pytest.importorskip("pyscf")
pytest.importorskip("basis_set_exchange")

from seamm_exec.testing import run_spec  # noqa: E402

pytestmark = pytest.mark.slow

SOURCE = Path(__file__).resolve().parent.parent

HI_STEP = textwrap.dedent('''
    """A step for testing: HI at the experimental r_e (1.6092 A)."""
    import seamm

    class HydrogenIodide(seamm.Node):
        def __init__(self, flowchart=None, extension=None):
            super().__init__(flowchart=flowchart, title="HI", extension=extension)

        @property
        def version(self):
            return "0.1"

        def description_text(self, P=None):
            return self.header + "\\n    HI at r = 1.6092 A."

        def run(self):
            next_node = super().run(None)
            db = self.get_variable("_system_db")
            system = db.create_system(name="HI")
            configuration = system.create_configuration(name="HI")
            configuration.atoms.append(
                x=[0.0, 0.0], y=[0.0, 0.0], z=[0.0, 1.6092], symbol=["H", "I"]
            )
            db.system = system
            return next_node

    class HydrogenIodideStep:
        my_description = {"description": "HI for tests", "group": "Building",
                          "name": "HydrogenIodide"}

        def __init__(self, flowchart=None, gui=None):
            pass

        def description(self):
            return HydrogenIodideStep.my_description

        def create_node(self, flowchart=None, **kwargs):
            return HydrogenIodide(flowchart=flowchart, **kwargs)

        def create_tk_node(self, canvas=None, **kwargs):
            raise NotImplementedError("no GUI for the test step")
''')

SPEC = """
title: HI spin-orbit correction
steps:
- HydrogenIodide: {}
- PySCF:
    steps:
    - SpinOrbit
"""


def install_hi_step(site):
    site.mkdir(parents=True, exist_ok=True)
    info = site / "seamm_test_hi-0.1.dist-info"
    info.mkdir(exist_ok=True)
    (site / "seamm_test_hi.py").write_text(HI_STEP)
    (info / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: seamm-test-hi\nVersion: 0.1\n"
    )
    (info / "entry_points.txt").write_text(
        "[org.molssi.seamm]\nHydrogenIodide = seamm_test_hi:HydrogenIodideStep\n\n"
        "[org.molssi.seamm.tk]\nHydrogenIodide = seamm_test_hi:HydrogenIodideStep\n"
    )
    return site


def test_hi_spin_orbit_correction(tmp_path):
    """HI: delta_SO 2.190 kJ/mol and the I splitting 7560 cm^-1, as in SEAMM 5877."""
    site = install_hi_step(tmp_path / "hi_site")
    job = run_spec(
        tmp_path / "run",
        SPEC,
        inis={"pyscf": {"code": sys.executable}},
        extra_path=[site],
        timeout=900,
    )
    text = "\n".join(p.read_text() for p in job.rglob("step.out"))
    text = re.sub(r"\s+", " ", text)
    m = re.search(r"E_SO\(J average\) = (-?[0-9.]+) kJ/mol", text)
    assert m, text[-3000:]
    assert float(m.group(1)) == pytest.approx(2.190, abs=0.005)
    m = re.search(r"I atom \(ECP28MDF, unc-cc-pVTZ-PP\): J splitting ([0-9.]+)", text)
    assert m and float(m.group(1)) == pytest.approx(7560.4, abs=0.5)
    assert "Known offset" in text
