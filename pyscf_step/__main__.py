# -*- coding: utf-8 -*-

"""Handle the installation of the PySCF step."""

from .installer import Installer


def run():
    """Create the conda environment `seamm-pyscf` with PySCF and write pyscf.ini."""
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
