***************
Getting Started
***************

Installation
============
The PySCF step is installed with the `SEAMM Manager`_. To add it, or bring it up to
date::

  seamm-manager install pyscf-step
  seamm-manager update pyscf-step

or use the Manager's window.

The calculations are done by PySCF_, which lives in its own conda environment,
``seamm-pyscf``. Installing the step creates it: the Manager runs the step's
installer, which

* installs PySCF 2.15 from PyPI (for the CPU; wheels exist for Linux and macOS);
* checks that PySCF imports and has the spin-orbit pseudopotential integrals the step
  uses;
* records the environment in ``pyscf.ini`` in the SEAMM installation's directory.

GPU4PySCF, for scalar HF and DFT on NVIDIA GPUs, is added only when asked for, with
``gpu = yes`` in ``pyscf.ini`` or ``--gpu``; the relativistic calculations always run on
the CPU. The installer can also be run by hand, e.g. to check a machine::

  pyscf-step-installer install
  pyscf-step-installer check
  pyscf-step-installer update

On a cluster, run it on a login node or where the jobs will run. A ``pyscf.ini`` with
``installation = local`` and ``code`` set to the full path of a Python that has PySCF
uses that environment instead, which the installer then leaves alone.

.. _SEAMM Manager: https://molssi-seamm.github.io/getting_started/installation/index.html
.. _PySCF: https://pyscf.org

A first calculation
===================
The Spin Orbit sub-step gives the molecular spin-orbit correction to the atomization
energy of a molecule containing heavy p-block atoms such as iodine. A flowchart for
hydrogen iodide::

  title: HI spin-orbit correction
  steps:
  - FromSMILESStep: {smiles string: I}
  - PySCF:
      steps:
      - SpinOrbit

Run it as a job. The step computes the iodine atom and the molecule, each with and
without spin-orbit, and prints the correction, about 2.2 kJ/mol for HI at its
experimental bond length. The correction is stored as the property
``spin-orbit correction#PySCF#HF/SO-ECP/unc-cc-pVTZ-PP``, so other steps and tables
can use it.

For what the numbers mean, the choices available and the method's validation, see the
:ref:`User Guide <user-guide>`.
