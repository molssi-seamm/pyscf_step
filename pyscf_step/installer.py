# -*- coding: utf-8 -*-

"""Installer for the PySCF plug-in.

PySCF runs in its own conda environment, ``seamm-pyscf`` by default, which this
installer creates and registers in ``pyscf.ini``. The environment holds Python 3.12
from conda-forge and PySCF from PyPI (wheels for Linux x86_64 and macOS arm64), for
the CPU. GPU4PySCF, for scalar HF/DFT on NVIDIA GPUs, is added only when
``gpu = yes`` is set in ``pyscf.ini`` (or ``--gpu`` is given); the relativistic
calculations always run on the CPU.
"""

import importlib.resources
import json
import logging
import os
import platform

import seamm_manager

logger = logging.getLogger(__name__)

#: The PySCF the step was validated with (SEAMM 5876/5877/5879, 2026-11-20). Pinned:
#: 2.15.1 has no macOS arm64 wheel.
PYSCF = "pyscf==2.15.0"

#: The optional GPU packages (CUDA 12 build, as the Interaction EDA step uses).
GPU_PACKAGES = ("gpu4pyscf-cuda12x>=1.9.1", "cupy-cuda12x[ctk]", "cutensor-cu12")

#: Imports PySCF and checks the pieces the step uses; no calculation.
PROBE = r"""
import json
out = {}
try:
    import pyscf
    from pyscf import gto, scf
    out["pyscf"] = pyscf.__version__
    out["ghf_soc"] = hasattr(scf.ghf.GHF, "with_soc")
    out["ecpso"] = "ECPso" in gto.moleintor._INTOR_FUNCTIONS
except Exception as e:
    out["import_error"] = f"{type(e).__name__}: {e}"
try:
    import gpu4pyscf
    out["gpu4pyscf"] = gpu4pyscf.__version__
except Exception:
    pass
print("PROBE" + json.dumps(out))
"""


def parse(output, tag):
    """The JSON after ``tag`` in a script's output, or None."""
    for line in (output or "").splitlines():
        if line.startswith(tag):
            try:
                return json.loads(line[len(tag) :])
            except ValueError:
                return None
    return None


def assess(probe):
    """(usable, [lines]) for a probe of the environment."""
    if not isinstance(probe, dict):
        return False, ["!   Could not run Python in the environment."]
    if "import_error" in probe:
        return False, [f"!   PySCF does not import: {probe['import_error']}"]
    lines = [f"    PySCF {probe['pyscf']}."]
    ok = True
    if not (probe.get("ghf_soc") and probe.get("ecpso")):
        ok = False
        lines.append(
            "!   This PySCF lacks spin-orbit ECP integrals in GHF; 2.15 is needed."
        )
    if "gpu4pyscf" in probe:
        lines.append(f"    GPU4PySCF {probe['gpu4pyscf']} is also installed.")
    return ok, lines


class Installer(seamm_manager.InstallerBase):
    """Handle the installation of PySCF for the PySCF step.

    #. ``install`` creates the conda environment ``seamm-pyscf`` (Python 3.12),
       installs PySCF in it (and GPU4PySCF if asked), checks it and writes
       ``pyscf.ini``.
    #. ``update`` upgrades PySCF (within the validated series) in that environment.
    #. ``check`` and ``show`` say what is installed and whether it works.

    An environment named in ``pyscf.ini`` that SEAMM did not create, or an
    ``installation = local`` Python, is never changed.
    """

    def __init__(self, logger=logger):
        super().__init__(logger=logger)
        self.environment = "seamm-pyscf"
        self.section = "pyscf-step"
        self.init_file_name = "pyscf.ini"
        self.executables = ["python"]
        self.resource_path = importlib.resources.files("pyscf_step") / "data"
        self.environment_file = self.resource_path / "seamm-pyscf.yml"

    def setup_parser(self):
        parser = super().setup_parser()
        for name in ("install", "update"):
            self.subparser[name].add_argument(
                "--gpu",
                action="store_true",
                help="Also install GPU4PySCF (Linux with an NVIDIA GPU only).",
            )
        return parser

    # ------------------------------------------------------------------
    def _packages(self):
        packages = [PYSCF]
        gpu = getattr(self.options, "gpu", False)
        if not gpu and self.exe_config.path is not None:
            value = str(self.exe_config.get_values("local").get("gpu", "no"))
            gpu = value.strip().lower() in ("yes", "true", "1")
        if gpu:
            if platform.system() != "Linux":
                print(
                    "    GPU4PySCF needs Linux and an NVIDIA GPU; it is not installed."
                )
            else:
                packages.extend(GPU_PACKAGES)
        return packages, gpu

    def _finish(self, environment):
        ok, lines = assess(parse(self.conda.run_python(environment, PROBE), "PROBE"))
        for line in lines:
            print(line)
        return ok

    def _register(self, environment, gpu):
        self.check_exe_configuration_file()
        self.exe_config.set_value("local", "installation", "conda")
        conda_exe = seamm_manager.find_conda()
        if conda_exe is not None:
            self.exe_config.set_value("local", "conda", conda_exe)
        self.exe_config.set_value("local", "conda-environment", environment)
        self.exe_config.set_value("local", "code", "python")
        self.exe_config.set_value("local", "gpu", "yes" if gpu else "no")
        for key in ("modules", "container", "platform"):
            self.exe_config.set_value("local", key, None)
        self.exe_config.save()

    def _exists(self, environment):
        if os.path.isabs(os.path.expanduser(environment)):
            return os.path.isdir(
                os.path.join(os.path.expanduser(environment), "conda-meta")
            )
        return self.conda.exists(environment)

    # ------------------------------------------------------------------
    def install(self):
        """Create the environment and install PySCF."""
        if self.shared_codes:
            self._use_shared_code()
            return
        if self.code_environments == "prefixed":
            from seamm_manager.naming import compute_tag
            from seamm_manager.policy import prefixed_environment

            self.environment = prefixed_environment(
                self.environment, compute_tag(self.root)
            )
        self.check_exe_configuration_file()
        data = self.exe_config.get_values("local")
        if data.get("installation") == "local":
            print(
                f"    pyscf.ini uses the local Python '{data.get('code')}', which is "
                "left as it is."
            )
            return
        packages, gpu = self._packages()
        environment = data.get("conda-environment") or self.environment
        if self._exists(environment):
            print(f"    Using the existing Conda environment '{environment}'.")
            if self._not_ours(environment):
                print(
                    f"    '{environment}' was not created by SEAMM, so nothing is "
                    "installed in it."
                )
                self._finish(environment)
                return
        else:
            print(f"    Installing Conda environment '{environment}' with PySCF.")
            self.conda.create_environment(
                self.environment_file, name=environment, conda_only=True
            )
        self.conda.pip_install(environment, packages)
        self._register(environment, gpu)
        if self._finish(environment):
            print("    Done!\n")
        else:
            print(
                "!   The environment was made but does not work yet; see above. "
                "'pyscf-step-installer check' checks it again.\n"
            )

    def update(self):
        """Upgrade PySCF in the environment."""
        if self.shared_codes:
            self._use_shared_code()
            return
        self.check_exe_configuration_file()
        data = self.exe_config.get_values("local")
        if data.get("installation") != "conda":
            print("!   PySCF was not installed by SEAMM, so it is not updated.")
            return
        environment = data.get("conda-environment") or self.environment
        if not self._exists(environment):
            print(f"!   The Conda environment '{environment}' does not exist; install.")
            return
        if self._not_ours(environment):
            print(f"    '{environment}' was not created by SEAMM; left as it is.")
            return
        packages, gpu = self._packages()
        print(f"    Updating PySCF in '{environment}'.")
        self.conda.pip_install(environment, packages, upgrade=True)
        self._register(environment, gpu)
        if self._finish(environment):
            print("    Done!\n")

    def check(self):
        """Check that pyscf.ini points at a working PySCF."""
        if not self.configuration.section_exists(self.section):
            self.check_configuration_file()
        self.check_exe_configuration_file()
        data = self.exe_config.get_values("local")
        if data.get("installation") == "conda":
            environment = data.get("conda-environment") or self.environment
            if not self._exists(environment):
                print(
                    f"!   The Conda environment '{environment}' in pyscf.ini does not "
                    "exist. Run 'pyscf-step-installer install'."
                )
                return False
            print(f"    The Conda environment is '{environment}'.")
            return self._finish(environment)
        print(f"    pyscf.ini runs '{data.get('code')}', which is not checked.")
        return True

    def exe_version(self, config):
        """The version of PySCF in its environment."""
        environment = config.get("conda-environment") or self.environment
        probe = parse(self.conda.run_python(environment, PROBE), "PROBE") or {}
        return "PySCF", probe.get("pyscf", "unknown")
