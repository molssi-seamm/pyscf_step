2026-11-20 Initial development: the spin-orbit correction
=========================================================

Goal of v1
----------
A step giving the molecular spin-orbit correction for thermochemistry of molecules with
iodine and other heavy p-block elements, at the level Bross et al. (J. Phys. Chem. A
127, 704 (2023)) use for the Active Thermochemical Tables: two-component SCF with
spin-orbit pseudopotentials minus scalar SCF. ORCA cannot do it: it has no spin-orbit
coupling with pseudopotentials and silently drops their spin-orbit part.

Decisions
---------
* **Spin-orbit pseudopotentials from scalar data.** The Basis Set Exchange and PySCF's
  library store the Stuttgart-Cologne MDF potentials as scalar terms only, but as
  j-weighted pairs (j = l - 1/2, l + 1/2, weights l/(2l+1), (l+1)/(2l+1)), so the
  spin-orbit operator U_SO l.s, U_SO = 2/(2l+1) [V(l, l+1/2) - V(l, l-1/2)], is -2c/l
  and +2c/(l+1) on the two primitives. This equals the published (Molpro-format) SO
  blocks to rounding; each element's data is checked against them
  (``data/so_ecp/stuttgart_so.json``).
* **Atoms by spherical average-of-configuration SCF only.** A first test with single
  determinants (SEAMM 5876) gave HI 2.80 kJ/mol against Bross's 2.4, and the same setup
  gave atomic splittings of 7868 and 7990 cm^-1 in two runs: broken-symmetry
  determinants relax differently with and without spin-orbit and are not reproducible.
  With aoc-SCF (5877) the I splitting is 7560 cm^-1 and HI 2.19 kJ/mol. PySCF has no
  aoc-SCF; ``engine/aoc.py`` implements it for GHF, checked against explicit
  determinant averages.
* **The correction includes the atomic J-average term.** The explicit 4s4p4d shells of a
  small-core pseudopotential carry a large second-order spin-orbit energy (~100 kJ/mol
  for I) that cancels between atom and molecule; alone it depends strongly on the basis
  set, in the correction it does not (HI 2.190 / 2.189 / 2.199 kJ/mol with
  unc-cc-pVTZ-PP / unc-aug-cc-pVTZ-PP / unc-cc-pVQZ-PP).
* **A known offset, not corrected:** 7-10 % below Bross et al. (CH3I 2.33 vs 2.6, I2
  7.57 vs 8.0 kJ/mol).
* **CPU only for the relativistic work.** GPU4PySCF 1.9.1's GHF has the spin-orbit
  pseudopotential terms but no spin-orbit X2C; for these molecule sizes the GPU gains
  nothing.

Tests
-----
The engine is checked against SEAMM 5877 (``tests/test_so_regression.py``) and the whole
sub-step by a flowchart run as a job (``tests/test_spin_orbit_run.py``); both need PySCF
and run in CI.
