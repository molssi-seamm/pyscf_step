===================
SEAMM PySCF Plug-in
===================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/pyscf_step
   :target: https://github.com/molssi-seamm/pyscf_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/pyscf_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/pyscf_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/pyscf_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/pyscf_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/pyscf_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/pyscf_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/pyscf_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/pyscf_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/pyscf_step.svg
   :target: https://pypi.python.org/pypi/pyscf_step
   :alt: PyPi VERSION

A SEAMM plug-in for PySCF

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/pyscf_step/index.html
* Code: https://github.com/molssi-seamm/pyscf_step

Features
--------

* The molecular spin-orbit correction to atomization energies (two-component minus
  scalar SCF with the spin-orbit part of energy-consistent pseudopotentials, the atoms
  by average-of-configuration SCF), for thermochemistry of molecules with I, Br and
  other heavy p-block elements.
* PySCF runs in its own environment, ``seamm-pyscf``, on the CPU, as SEAMM tasks.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
