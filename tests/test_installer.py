#!/usr/bin/env python
"""The installer's resource files ship with the package."""

import importlib.resources


def test_installer_resources():
    data = importlib.resources.files("pyscf_step") / "data"
    for name in ("configuration.txt", "pyscf.ini", "seamm-pyscf.yml"):
        assert (data / name).is_file(), name
    assert (data / "configuration.txt").read_text().startswith("[pyscf-step]")


def test_probe_assessment():
    from pyscf_step.installer import assess

    ok, _ = assess({"pyscf": "2.15.0", "ghf_soc": True, "ecpso": True})
    assert ok
    ok, lines = assess({"pyscf": "2.14.0", "ghf_soc": True, "ecpso": False})
    assert not ok and "2.15" in lines[-1]
    ok, _ = assess({"import_error": "ModuleNotFoundError: pyscf"})
    assert not ok
