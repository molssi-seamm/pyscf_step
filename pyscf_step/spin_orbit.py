# -*- coding: utf-8 -*-

"""Non-graphical part of the Spin Orbit step in a PySCF flowchart"""

import json
import logging
from pathlib import Path
import importlib.resources
import pprint  # noqa: F401

import pyscf_step  # noqa: E999
from pyscf_step import basis as so_basis
import molsystem
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("PySCF")

# Add this module's properties to the standard properties
path = importlib.resources.files("pyscf_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class SpinOrbit(pyscf_step.Substep):
    """
    The non-graphical part of a Spin Orbit step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : SpinOrbitParameters
        The control parameters for Spin Orbit.

    See Also
    --------
    TkSpinOrbit,
    SpinOrbit, SpinOrbitParameters
    """

    def __init__(
        self, flowchart=None, title="Spin Orbit", extension=None, logger=logger
    ):
        """A substep for Spin Orbit in a subflowchart for PySCF.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Spin Orbit {self}")

        super().__init__(
            flowchart=flowchart,
            title="Spin Orbit",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self._calculation = "Spin Orbit"
        self._model = None
        self._metadata = pyscf_step.metadata
        self.parameters = pyscf_step.SpinOrbitParameters()

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def version(self):
        """The semantic version of this module."""
        return pyscf_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return pyscf_step.__git_revision__

    def description_text(self, P=None):
        """A text description of what this step will do."""
        if not P:
            P = self.parameters.values_to_dict()
        unc = "uncontracted " if P["uncontract"] else ""
        text = (
            "The molecular spin-orbit correction to the atomization energy: "
            f"{P['method']} with and without the spin-orbit part of the "
            f"pseudopotentials ({P['route']}), {unc}{P['basis']} "
            f"({P['other basis']} for other elements). The atoms with spin-orbit "
            "pseudopotentials are done by spherical average-of-configuration SCF, "
            "so that the correction is relative to their J-averaged energies."
        )
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run the Spin Orbit step: the atoms, then the molecule, then delta_SO."""
        import seamm_exec

        next_node = super().run(printer)
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        _, configuration = self.get_system_configuration(None)
        if configuration.periodicity != 0:
            raise RuntimeError("The spin-orbit correction is for molecules only.")
        symbols = list(configuration.atoms.symbols)
        xyz = configuration.atoms.get_coordinates(fractionals=False)
        charge = int(configuration.charge)
        multiplicity = int(configuration.spin_multiplicity)
        if multiplicity != 1:
            raise RuntimeError(
                "The spin-orbit correction is implemented for closed-shell molecules "
                f"only (this one has multiplicity {multiplicity})."
            )

        basis, ecp, info = so_basis.basis_and_ecp(
            symbols, P["basis"], P["other basis"], uncontract=bool(P["uncontract"])
        )
        model = (
            f"{P['method']}/{P['route']}/"
            + ("unc-" if P["uncontract"] else "")
            + P["basis"]
        )
        self._model = model
        so_elements = sorted(ecp)
        indent = self.indent + 4 * " "
        if not so_elements:
            printer.normal(
                __(
                    "No element has a spin-orbit pseudopotential, so the correction "
                    "is zero.",
                    indent=indent,
                )
            )
            self.store_results(
                configuration=configuration,
                data={"spin-orbit correction": 0.0},
                create_tables=self.parameters["results"].value != {},
            )
            return next_node
        for element in so_elements:
            if element not in so_basis.P_SHELL:
                raise NotImplementedError(
                    f"The {element} atom's J states need a spinor CI (open p2, p3 or "
                    "p4 shell), which is not available yet."
                )

        # The tasks: one per element with an SO-ECP, and the molecule
        files = {
            name: (
                importlib.resources.files("pyscf_step") / "engine" / name
            ).read_text()
            for name in ("so_driver.py", "aoc.py", "so_correction.py")
        }
        task_set = seamm_exec.TaskSet(self, directory=directory)
        inputs = {}
        for element in so_elements:
            inputs[f"atom_{element}"] = {
                "task": "atom",
                "atoms": [[element, [0.0, 0.0, 0.0]]],
                "basis": {element: basis[element]},
                "ecp": {element: ecp[element]},
                "spin": 1,
                "n_p": so_basis.P_SHELL[element],
            }
        inputs["molecule"] = {
            "task": "molecule",
            "atoms": [[s, [float(x) for x in r]] for s, r in zip(symbols, xyz)],
            "basis": basis,
            "ecp": ecp,
            "charge": charge,
            "spin": 0,
        }
        for key, data in inputs.items():
            nao = sum(n_functions(basis[a[0]]) for a in data["atoms"])
            cores = 4 if data["task"] == "atom" else 8
            memory = max(4e9, 3.0 * nao**4)
            data["memory_mb"] = int(0.8 * memory / 1e6)
            data["log"] = "pyscf.log"
            task_files = dict(files)
            task_files["input.json"] = json.dumps(data, indent=1)
            task_set.add(
                seamm_exec.Task(
                    key=key,
                    program="pyscf",
                    cmd=["{code}", "so_driver.py", "input.json", "-o", "results.json"],
                    files=task_files,
                    return_files=["results.json", "pyscf.log"],
                    resources=seamm_exec.Resources(
                        ntasks=1, cpus_per_task=cores, mem_per_cpu=int(memory / cores)
                    ),
                    shell=True,
                )
            )

        results = {}
        for result in task_set.run():
            record = None
            if result.directory is not None:
                path = Path(result.directory) / "results.json"
                if path.exists():
                    record = json.loads(path.read_text())
            if record is None or "error" in record:
                reason = (record or {}).get("error") or result.reason or "no result"
                raise RuntimeError(f"PySCF task {result.key} failed: {reason}")
            results[result.key] = record

        mol = results["molecule"]
        atoms = {e: results[f"atom_{e}"] for e in so_elements}
        counts = {e: symbols.count(e) for e in so_elements}
        e_so_atoms = sum(counts[e] * atoms[e]["E_SO_Javg_kJ"] for e in so_elements)
        delta = -mol["E_SO_kJ"] + e_so_atoms
        data = {
            "spin-orbit correction": delta,
            "E_SO molecule": mol["E_SO_kJ"],
            "E_SO atoms": e_so_atoms,
            "scalar energy": mol["E_scalar"],
            "two-component energy": mol["E_2c"],
            "atomic spin-orbit": {
                e: {
                    "splitting_cm": atoms[e]["splitting_cm"],
                    "E_SO_ground_kJ": atoms[e]["E_SO_ground_kJ"],
                    "E_SO_Javg_kJ": atoms[e]["E_SO_Javg_kJ"],
                }
                for e in so_elements
            },
        }
        self.analyze(data=data, info=info, counts=counts, mol=mol)
        self.store_results(
            configuration=configuration,
            data=data,
            create_tables=self.parameters["results"].value != {},
        )
        self.references.cite(
            raw=BROSS_2023,
            alias="Bross2023",
            module="pyscf_step",
            level=1,
            note="The spin-orbit correction as two-component minus scalar SCF.",
        )
        self.references.cite(
            raw=PYSCF_2020,
            alias="pyscf2020",
            module="pyscf_step",
            level=1,
            note="The PySCF program package.",
        )
        for element in so_elements:
            for alias in ECP_PAPERS.get(element, ()):
                self.references.cite(
                    raw=ECP_BIBTEX[alias],
                    alias=alias,
                    module="pyscf_step",
                    level=2,
                    note=f"The {info[element]['ecp']} pseudopotential for {element}.",
                )
        return next_node

    def analyze(self, indent="", data=None, info=None, counts=None, mol=None, **kwargs):
        """Print the results."""
        if data is None:
            return
        indent = self.indent + 4 * " "
        lines = []
        for element, terms in data["atomic spin-orbit"].items():
            lines.append(
                f"{element} atom ({info[element]['ecp']}, {info[element]['basis']}): "
                f"J splitting {terms['splitting_cm']:.1f} cm^-1, "
                f"E_SO(ground J) {terms['E_SO_ground_kJ']:.3f} kJ/mol, "
                f"E_SO(J average) {terms['E_SO_Javg_kJ']:.3f} kJ/mol."
            )
        lines.append(
            f"Molecule: E(two-component) - E(scalar) = {data['E_SO molecule']:.3f} "
            f"kJ/mol; the scalar-limit check (GHF without spin-orbit - RHF) is "
            f"{mol.get('GHF_noSOC_minus_RHF', float('nan')):.1e} E_h."
        )
        lines.append(
            "The molecular spin-orbit correction to the atomization energy, "
            "delta_SO = -E_SO(molecule) + sum of the atoms' E_SO(J average) = "
            f"{data['spin-orbit correction']:.3f} kJ/mol."
        )
        lines.append(
            "Known offset: with uncontracted cc-pVTZ-PP and experimental geometries "
            "this method gives values 7-10 % below Bross et al.'s DIRAC values "
            "(HI 2.19 vs 2.4, CH3I 2.33 vs 2.6, I2 7.57 vs 8.0 kJ/mol), likely from "
            "the geometries and DIRAC's implementation. It is reported, not corrected."
        )
        for line in lines:
            printer.normal(__(line, indent=indent))
            printer.normal("")


def n_functions(shells):
    """The number of spherical functions in a PySCF basis list."""
    return sum((2 * shell[0] + 1) * (len(shell[1]) - 1) for shell in shells)


BROSS_2023 = """@article{Bross2023,
  title = {Active Thermochemical Tables: Enthalpies of Formation of Bromo- and
    Iodo-Methanes, Ethenes and Ethynes},
  author = {Bross, David H. and Bacskay, George B. and Peterson, Kirk A. and
    Ruscic, Branko},
  journal = {J. Phys. Chem. A},
  volume = {127},
  pages = {704--723},
  year = {2023},
  doi = {10.1021/acs.jpca.2c07897}
}"""

PYSCF_2020 = """@article{pyscf2020,
  title = {Recent developments in the {PySCF} program package},
  author = {Sun, Qiming and others},
  journal = {J. Chem. Phys.},
  volume = {153},
  pages = {024109},
  year = {2020},
  doi = {10.1063/5.0006074}
}"""

#: The pseudopotential papers, as the Basis Set Exchange cites them.
ECP_BIBTEX = {
    "peterson2003b": """@article{peterson2003b,
  author = {Peterson, Kirk A. and Figgen, Detlev and Goll, Erich and Stoll, Hermann
    and Dolg, Michael},
  title = {Systematically convergent basis sets with relativistic pseudopotentials.
    {II}. Small-core pseudopotentials and correlation consistent basis sets for the
    post-d group 16-18 elements},
  journal = {J. Chem. Phys.},
  volume = {119},
  pages = {11113--11123},
  year = {2003},
  doi = {10.1063/1.1622924}
}""",
    "peterson2006a": """@article{peterson2006a,
  author = {Peterson, Kirk A. and Shepler, Benjamin C. and Figgen, Detlev and
    Stoll, Hermann},
  title = {On the Spectroscopic and Thermochemical Properties of {ClO}, {BrO}, {IO},
    and Their Anions},
  journal = {J. Phys. Chem. A},
  volume = {110},
  pages = {13877--13883},
  year = {2006},
  doi = {10.1021/jp065887l}
}""",
}

#: The papers for each element's pseudopotential.
ECP_PAPERS = {"Br": ("peterson2003b",), "I": ("peterson2003b", "peterson2006a")}
