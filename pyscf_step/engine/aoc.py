# -*- coding: utf-8 -*-

"""Average-of-configuration SCF (aoc-SCF) for atoms, in PySCF's GHF.

PySCF has no average-of-configuration SCF. It is needed for the spin-orbit
correction: single determinants of open-shell atoms break the spherical symmetry
and relax differently with and without spin-orbit, and the results are not even
reproducible (phase 0, SEAMM 5876/5877). The spherical aoc treatment, as Bross,
Bacskay, Peterson and Ruscic (J. Phys. Chem. A 127, 704 (2023)) use with DIRAC
(Visser et al., J. Chem. Phys. 96, 2910 (1992)), is the only atomic treatment here.

Validated in SEAMM 5877 (ECP28MDF with SO, uncontracted cc-pVTZ-PP): the I atom 2P
splitting is 7560.4 cm-1 (7603 measured; DIRAC aoc-SCF 7587).
"""

import itertools

import numpy as np

HARTREE_TO_KJ = 2625.499639
HARTREE_TO_CM = 219474.6313632


def _pyscf():
    from pyscf import lib, scf

    return lib, scf


def make_aoc_class():
    """The AOCGHF class, built on PySCF's GHF (imported only when needed)."""
    lib, scf = _pyscf()

    class AOCGHF(scf.ghf.GHF):
        """Average-of-configuration SCF for one open shell, in GHF.

        Closed orbitals have occupation 1; the m open orbitals hold n electrons,
        averaged over all C(m, n) determinants. With f = n/m, a = (n-1)/(m-1), unit
        projectors D_c, D_o, D1 = D_c + f D_o and G = J - K,

            E = Tr h D1 + 1/2 Tr D1 G(D1) + 1/2 f (a - f) Tr D_o G(D_o),

        i.e. the open-open pairs carry f a instead of fractional occupation's f^2.
        The orbitals are optimised with an effective Fock in the MO basis whose
        off-diagonal blocks are the energy gradient: closed-virtual F_c,
        open-virtual F_o, closed-open F_c - f (a - f)/(1 - f) G_o, with
        F_c = h + G(D1) and F_o = F_c + (a - f) G(D_o). The closed and open sets are
        followed by maximum overlap with the previous iteration's.
        """

        _keys = scf.ghf.GHF._keys | {
            "aoc_f",
            "aoc_a",
            "aoc_m",
            "aoc_n",
            "_ref_c",
            "_ref_o",
        }

        def __init__(self, mol, with_soc, c_closed, c_open, n_open):
            scf.ghf.GHF.__init__(self, mol)
            m = c_open.shape[1]
            if not 0 < n_open < m:
                raise ValueError(
                    f"an open shell needs 0 < n < m, not n={n_open}, m={m}"
                )
            self.with_soc = with_soc
            self._ref_c = c_closed
            self._ref_o = c_open
            self.aoc_m = m
            self.aoc_n = n_open
            self.aoc_f = n_open / m
            self.aoc_a = (n_open - 1) / (m - 1)

        def get_occ(self, mo_energy=None, mo_coeff=None):
            if mo_coeff is None:
                raise RuntimeError("aoc-SCF needs orbitals")
            s = self.get_ovlp()
            nc = self._ref_c.shape[1]
            ov = self._ref_c.conj().T @ s @ mo_coeff
            pc = np.einsum("ij,ij->j", ov.conj(), ov).real
            ic = np.argsort(-pc)[:nc]
            rest = np.setdiff1d(np.arange(mo_coeff.shape[1]), ic)
            ov = self._ref_o.conj().T @ s @ mo_coeff[:, rest]
            po = np.einsum("ij,ij->j", ov.conj(), ov).real
            io = rest[np.argsort(-po)[: self.aoc_m]]
            occ = np.zeros(mo_coeff.shape[1])
            occ[ic] = 1.0
            occ[io] = self.aoc_f
            self._ref_c = mo_coeff[:, ic]
            self._ref_o = mo_coeff[:, io]
            return occ

        def make_rdm1(self, mo_coeff=None, mo_occ=None, **kwargs):
            if mo_coeff is None:
                mo_coeff = self.mo_coeff
            if mo_occ is None:
                mo_occ = self.mo_occ
            dm = (mo_coeff * mo_occ) @ mo_coeff.conj().T
            return lib.tag_array(dm, mo_coeff=mo_coeff, mo_occ=mo_occ)

        def _parts(self, dm):
            c_all = getattr(dm, "mo_coeff", None)
            occ = getattr(dm, "mo_occ", None)
            if c_all is None:
                return None
            o = (occ > 1e-12) & (occ < 1 - 1e-12)
            c = occ >= 1 - 1e-12
            v = occ <= 1e-12
            c_o = c_all[:, o]
            d_o = c_o @ c_o.conj().T
            g_o = scf.ghf.GHF.get_veff(self, self.mol, d_o)
            return c_all, c, o, v, d_o, g_o

        def get_fock(
            self,
            h1e=None,
            s1e=None,
            vhf=None,
            dm=None,
            cycle=-1,
            diis=None,
            diis_start_cycle=None,
            level_shift_factor=None,
            damp_factor=None,
            fock_last=None,
            **kwargs,
        ):
            if h1e is None:
                h1e = self.get_hcore()
            if s1e is None:
                s1e = self.get_ovlp()
            if dm is None:
                dm = self.make_rdm1()
            if vhf is None:
                vhf = self.get_veff(self.mol, dm)
            f_c = h1e + vhf
            parts = self._parts(dm)
            if parts is None:
                fock = f_c
            else:
                c_all, c, o, v, d_o, g_o = parts
                f, a = self.aoc_f, self.aoc_a
                fc = c_all.conj().T @ f_c @ c_all
                go = c_all.conj().T @ g_o @ c_all
                fo = fc + (a - f) * go
                k = f * (a - f) / (1 - f)
                fm = fc.copy()
                fm[np.ix_(o, o)] = fo[np.ix_(o, o)]
                fm[np.ix_(o, v)] = fo[np.ix_(o, v)]
                fm[np.ix_(v, o)] = fo[np.ix_(v, o)]
                fm[np.ix_(c, o)] = fc[np.ix_(c, o)] - k * go[np.ix_(c, o)]
                fm[np.ix_(o, c)] = fc[np.ix_(o, c)] - k * go[np.ix_(o, c)]
                sc = s1e @ c_all
                fock = sc @ fm @ sc.conj().T
            start = (
                self.diis_start_cycle if diis_start_cycle is None else diis_start_cycle
            )
            if diis is not None and cycle >= start:
                fock = diis.update(s1e, dm, fock, self, h1e, vhf)
            return fock

        def energy_elec(self, dm=None, h1e=None, vhf=None):
            if dm is None:
                dm = self.make_rdm1()
            if h1e is None:
                h1e = self.get_hcore()
            if vhf is None:
                vhf = self.get_veff(self.mol, dm)
            e1 = np.einsum("ij,ji->", h1e, dm).real
            e2 = 0.5 * np.einsum("ij,ji->", vhf, dm).real
            parts = self._parts(dm)
            if parts is not None:
                d_o, g_o = parts[4], parts[5]
                e2 += (
                    0.5
                    * self.aoc_f
                    * (self.aoc_a - self.aoc_f)
                    * np.einsum("ij,ji->", g_o, d_o).real
                )
            self.scf_summary["e1"] = e1
            self.scf_summary["e2"] = e2
            return e1 + e2, e2

        def get_grad(self, mo_coeff, mo_occ, fock=None):
            if fock is None:
                fock = self.get_fock(dm=self.make_rdm1(mo_coeff, mo_occ))
            fm = mo_coeff.conj().T @ fock @ mo_coeff
            o = (mo_occ > 1e-12) & (mo_occ < 1 - 1e-12)
            c = mo_occ >= 1 - 1e-12
            v = mo_occ <= 1e-12
            return np.concatenate(
                [
                    fm[np.ix_(c, v)].ravel(),
                    fm[np.ix_(o, v)].ravel(),
                    fm[np.ix_(c, o)].ravel(),
                ]
            )

    return AOCGHF


def determinant_average(mol, with_soc, mf):
    """The average of the configuration's determinant energies at mf's orbitals.

    A check of the aoc energy expression: equal to mf.e_tot to rounding.
    """
    lib, scf = _pyscf()
    plain = scf.GHF(mol)
    plain.with_soc = with_soc
    c_all, occ = mf.mo_coeff, mf.mo_occ
    c_idx = np.where(occ >= 1 - 1e-12)[0]
    o_idx = np.where((occ > 1e-12) & (occ < 1 - 1e-12))[0]
    d_c = c_all[:, c_idx] @ c_all[:, c_idx].conj().T
    energies = []
    for combo in itertools.combinations(o_idx, mf.aoc_n):
        cc = c_all[:, list(combo)]
        energies.append(plain.energy_tot(d_c + cc @ cc.conj().T))
    return float(np.mean(energies)), len(energies)


def _fractional_occ(n_full, n_open, n_el_open):
    def get_occ(mo_energy=None, mo_coeff=None):
        idx = np.argsort(mo_energy)
        occ = np.zeros(len(mo_energy))
        occ[idx[:n_full]] = 1.0
        occ[idx[n_full : n_full + n_open]] = n_el_open / n_open
        return occ

    return get_occ


def _ghf_dm(dm):
    dm = np.asarray(dm)
    da, db = (dm / 2, dm / 2) if dm.ndim == 2 else dm
    n = da.shape[0]
    out = np.zeros((2 * n, 2 * n), dtype=complex)
    out[:n, :n] = da
    out[n:, n:] = db
    return out


def _converge(mf, dm0, conv_tol=1e-10, max_cycle=300):
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    e = mf.kernel(dm0)
    if not mf.converged:
        raise RuntimeError(f"{type(mf).__name__} did not converge")
    return e


def p_shell_atom(mol, n_p, check=True, log=None):
    """Scalar and two-component aoc energies of an atom with an open p shell.

    For p1 and p5 the jj configurations are exact J states at SCF level: p1 gives
    J = 1/2 (p1/2^1) and J = 3/2 (p3/2^1); p5 gives J = 3/2 (p1/2^2 p3/2^3) and
    J = 1/2 (p1/2^1 p3/2^4). Each state is the aoc over its jj determinants with its
    own orbitals; the scalar 2P is the aoc over the six determinants of p^n. Other
    p shells (p2, p3, p4) need a small spinor CI and are refused.

    The p shell is taken as the six highest occupied spin-orbitals of a spherical,
    fractionally occupied calculation, which must split 2 + 4 (p1/2, p3/2) with
    spin-orbit and be degenerate without it.

    Returns a dictionary of energies (Eh), the splitting (cm-1) and the spin-orbit
    terms (kJ/mol): E_SO(ground) = E(ground J) - E(2P) and
    E_SO(J-avg) = E(J-avg) - E(2P), with E(J-avg) the (2J+1)-weighted average.
    """
    if n_p not in (1, 5):
        raise NotImplementedError(
            f"p{n_p} atoms need a spinor CI for their J states; only p1 and p5 are "
            "available"
        )
    lib, scf = _pyscf()
    AOCGHF = make_aoc_class()
    say = log or (lambda *a: None)
    n_core = mol.nelectron - n_p
    out = {}
    for with_soc in (False, True):
        uhf = scf.UHF(mol)
        _converge(uhf, None)
        sph = scf.GHF(mol)
        sph.with_soc = with_soc
        sph.get_occ = _fractional_occ(n_core, 6, n_p)
        _converge(sph, _ghf_dm(uhf.make_rdm1()))
        order = np.argsort(sph.mo_energy)
        core, p = order[:n_core], order[n_core : n_core + 6]
        e_p = sph.mo_energy[p]
        c_sph = sph.mo_coeff
        if with_soc:
            if not (
                abs(e_p[1] - e_p[0]) < 1e-5
                and np.ptp(e_p[2:]) < 1e-5
                and e_p[2] - e_p[1] > 1e-3
            ):
                raise RuntimeError(f"the p spinors do not split 2 + 4: {e_p}")
            p12, p32 = p[:2], p[2:]
            if n_p == 5:
                configs = {
                    "J=3/2": (np.r_[core, p12], p32, 3),
                    "J=1/2": (np.r_[core, p32], p12, 1),
                }
            else:
                configs = {
                    "J=1/2": (core, p12, 1),
                    "J=3/2": (core, p32, 1),
                }
        else:
            if np.ptp(e_p) > 1e-5:
                raise RuntimeError(
                    f"the scalar p spin-orbitals are not degenerate: {e_p}"
                )
            configs = {"2P": (core, p, n_p)}
        for name, (closed, opened, n) in configs.items():
            mf = AOCGHF(mol, with_soc, c_sph[:, closed], c_sph[:, opened], n)
            occ0 = np.zeros(c_sph.shape[1])
            occ0[closed] = 1.0
            occ0[opened] = n / len(opened)
            e = _converge(mf, mf.make_rdm1(c_sph, occ0))
            entry = {"E": float(e)}
            if check:
                e_avg, ndet = determinant_average(mol, with_soc, mf)
                entry["determinant_average_check"] = abs(e_avg - e)
                entry["n_determinants"] = ndet
                if abs(e_avg - e) > 1e-8:
                    raise RuntimeError(
                        f"aoc energy check failed for {name}: {abs(e_avg - e):.2e}"
                    )
            say(f"aoc {name}: E = {e:.10f}")
            out[name] = entry
    e_sc, e12, e32 = out["2P"]["E"], out["J=1/2"]["E"], out["J=3/2"]["E"]
    e_javg = (2 * e12 + 4 * e32) / 6
    ground = "J=3/2" if n_p == 5 else "J=1/2"
    e_ground = out[ground]["E"]
    splitting = abs(e12 - e32) * HARTREE_TO_CM
    return {
        "states": out,
        "ground": ground,
        "E_scalar": e_sc,
        "E_Javg": e_javg,
        "splitting_cm": splitting,
        "E_SO_ground_kJ": (e_ground - e_sc) * HARTREE_TO_KJ,
        "E_SO_Javg_kJ": (e_javg - e_sc) * HARTREE_TO_KJ,
    }
