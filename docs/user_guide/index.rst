.. _user-guide:

**********
User Guide
**********

The PySCF step runs PySCF_ calculations in a SEAMM flowchart. Like the other code
steps it holds sub-steps; this version has one, **Spin Orbit**, for the spin-orbit
correction in thermochemistry. Each calculation is a task run in the ``seamm-pyscf``
environment, on the job's machine or, with a task queue, on a cluster, like any other
code.

.. _PySCF: https://pyscf.org

Spin Orbit
==========

What it computes
----------------
Composite thermochemistry schemes treat a molecule as free of spin-orbit coupling and
subtract the experimental spin-orbit lowering of the atoms (the J-averaged value, 30.3
kJ/mol for iodine). A closed-shell molecule is not free of it: spin-orbit coupling mixes
excited states into its ground state and lowers its energy. For molecules with iodine
this second-order term is a few kJ/mol, and leaving it out makes their heats of
formation too high by about that much (2.7-3.1 kJ/mol for HI and CH\ :sub:`3`\ I with
a CCSD(T) chain). It is larger for the heavier elements.

The step computes it as Bross, Bacskay, Peterson and Ruscic do for the Active
Thermochemical Tables [Bross2023]_: two-component SCF with the spin-orbit part of the
pseudopotentials, minus the same scalar calculation. For a molecule M,

.. math::

   E_\mathrm{SO}(\mathrm{M}) = E_\mathrm{2c}(\mathrm{M}) - E_\mathrm{scalar}(\mathrm{M})

   \delta_\mathrm{SO} = -E_\mathrm{SO}(\mathrm{M})
       + \sum_\mathrm{atoms} E_\mathrm{SO}(\mathrm{atom}, J\text{-average})

:math:`\delta_\mathrm{SO}` is the correction to the atomization energy that goes with
the experimental atomic J-averages; it is positive (the molecule is lowered) and is
added to the atomization energy, i.e. subtracted from the heat of formation. The atomic
term removes the large spin-orbit energy of the core electrons a small-core
pseudopotential leaves explicit (about 100 kJ/mol for iodine's 4s4p4d), which cancels
between atoms and molecule. Bross et al. found that the SCF level fits the ATcT data
best (RMS 0.79 kJ/mol over 11 molecules, against 3.6 without the term).

How
---
* **Pseudopotentials.** The basis sets and their energy-consistent (Stuttgart-Cologne)
  pseudopotentials come from the Basis Set Exchange. The libraries store only the scalar
  potentials; the spin-orbit part is rebuilt from their j-dependent terms and checked
  against the published spin-orbit parameters shipped with the step (iodine ECP28MDF
  and bromine ECP10MDF so far). An element with a pseudopotential but no such data is
  refused rather than run without spin-orbit.
* **The molecule.** Scalar RHF, and two-component GHF with the spin-orbit
  pseudopotential, at the same geometry. A check that GHF without spin-orbit equals RHF
  is printed.
* **The atoms.** Each element with a spin-orbit pseudopotential is computed by spherical
  average-of-configuration SCF, scalar and two-component alike: the scalar :sup:`2`\ P
  level, and the two J levels (for p\ :sup:`5`, p\ :sub:`1/2`\ :sup:`2`\
  p\ :sub:`3/2`\ :sup:`3` is J = 3/2 and p\ :sub:`1/2`\ :sup:`1`\ p\ :sub:`3/2`\ :sup:`4`
  is J = 1/2), each with its own orbitals. Single determinants are not used: they break
  the atom's symmetry, relax differently with and without spin-orbit and are not even
  reproducible. Atoms with p\ :sup:`1` or p\ :sup:`5` valence shells are available (the
  groups of B and F); those with p\ :sup:`2`, p\ :sup:`3` or p\ :sup:`4` (Sn, Sb, Te,
  Pb, Bi, ...) need a small spinor CI and are refused for now.

Parameters
----------
Relativistic route
    ``SO-ECP``, the spin-orbit pseudopotential route above.
Method
    ``HF``, the level Bross et al. recommend.
Basis set
    The basis set for the elements it covers; its pseudopotentials must have published
    spin-orbit data. The default is ``cc-pVTZ-PP``.
Basis set for other elements
    For the lighter elements, all-electron and without spin-orbit; default ``cc-pVTZ``.
Uncontract the basis sets
    ``yes`` (the default) uses them fully uncontracted, as Bross and Peterson recommend.

Results
-------
.. list-table::
   :header-rows: 1

   * - Result
     - Units
     - Meaning
   * - spin-orbit correction
     - kJ/mol
     - :math:`\delta_\mathrm{SO}`; also stored as the property
       ``spin-orbit correction#PySCF#{model}``, e.g. model
       ``HF/SO-ECP/unc-cc-pVTZ-PP``
   * - E_SO molecule
     - kJ/mol
     - :math:`E_\mathrm{2c} - E_\mathrm{scalar}` of the molecule
   * - E_SO atoms
     - kJ/mol
     - The sum of the atoms' :math:`E_\mathrm{SO}(J\text{-average})`
   * - scalar energy, two-component energy
     - E\ :sub:`h`
     - The molecule's two SCF energies
   * - atomic spin-orbit
     - cm\ :sup:`-1`, kJ/mol
     - Per element: the computed J splitting and :math:`E_\mathrm{SO}` of the ground J
       level and of the J average

Only closed-shell molecules are handled. The geometry is the current configuration's:
use the same one as for the rest of the thermochemistry.

Validation and the known offset
-------------------------------
With the defaults and experimental bond lengths:

.. list-table::
   :header-rows: 1

   * - Species
     - This step
     - Bross et al. (SCF)
   * - I atom :sup:`2`\ P splitting
     - 7560 cm\ :sup:`-1`
     - 7587 cm\ :sup:`-1` (experiment 7603)
   * - HI
     - 2.19 kJ/mol
     - 2.4
   * - CH\ :sub:`3`\ I
     - 2.33 kJ/mol
     - 2.6
   * - I\ :sub:`2`
     - 7.57 kJ/mol
     - 8.0

The basis set is not the cause of the differences: HI gives 2.190, 2.189 and 2.199
kJ/mol with uncontracted cc-pVTZ-PP, aug-cc-pVTZ-PP and cc-pVQZ-PP. Our values are
consistently 7-10 % below Bross et al.'s, likely from their CCSD(T)-F12b geometries and
the details of DIRAC's average-of-configuration and two-component pseudopotential
implementations. The offset is reported with every result and not corrected. ORCA's
all-electron CASSCF/NEVPT2 spin-orbit route gives 2.25-2.29 kJ/mol for HI.

References
----------
.. [Bross2023] D. H. Bross, G. B. Bacskay, K. A. Peterson and B. Ruscic, "Active
   Thermochemical Tables: Enthalpies of Formation of Bromo- and Iodo-Methanes, Ethenes
   and Ethynes", *J. Phys. Chem. A* **127**, 704-723 (2023),
   doi:10.1021/acs.jpca.2c07897.

The pseudopotentials: K. A. Peterson, D. Figgen, E. Goll, H. Stoll and M. Dolg, *J.
Chem. Phys.* **119**, 11113 (2003) (Br); K. A. Peterson, B. C. Shepler, D. Figgen and
H. Stoll, *J. Phys. Chem. A* **110**, 13877 (2006) (I). PySCF: Q. Sun et al., *J. Chem.
Phys.* **153**, 024109 (2020).

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
