=======
History
=======

2026.11.20 -- First release: the spin-orbit correction for thermochemistry
    * The Spin Orbit sub-step gives the molecular (second-order) spin-orbit correction
      to the atomization energy: two-component SCF with the spin-orbit part of the
      energy-consistent pseudopotentials minus the same scalar calculation, with the
      atoms by spherical average-of-configuration SCF, as Bross et al. do for the
      Active Thermochemical Tables. It is stored as the property "spin-orbit
      correction#PySCF#{model}".
    * The spin-orbit pseudopotentials are rebuilt from the Basis Set Exchange's data
      and checked against the published parameters, for iodine and bromine so far;
      other elements with pseudopotentials are refused rather than run without
      spin-orbit. Atoms with p1 or p5 valence shells are available.
    * Validated against Bross et al. (2023): the I atom splitting is 7560 cm^-1
      (7587 in their work; 7603 measured) and HI's correction 2.19 kJ/mol (2.4); our
      values are 7-10 % below theirs, which is reported with each result.
    * PySCF runs in its own conda environment, seamm-pyscf, which the installer
      creates (PySCF 2.15.0, CPU; GPU4PySCF only on request).
