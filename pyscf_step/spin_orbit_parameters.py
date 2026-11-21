# -*- coding: utf-8 -*-
"""
Control parameters for the Spin Orbit step in a SEAMM flowchart
"""

import logging
import seamm
import pprint  # noqa: F401

logger = logging.getLogger(__name__)


class SpinOrbitParameters(seamm.Parameters):  # noqa: E999
    """The control parameters for the Spin Orbit sub-step.

    See Also
    --------
    SpinOrbit, TkSpinOrbit
    """

    parameters = {
        "route": {
            "default": "SO-ECP",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("SO-ECP",),
            "format_string": "s",
            "description": "Relativistic route:",
            "help_text": (
                "SO-ECP: two-component SCF with the spin-orbit part of the "
                "energy-consistent pseudopotentials, minus the same scalar "
                "calculation (Bross et al., J. Phys. Chem. A 127, 704 (2023))."
            ),
        },
        "method": {
            "default": "HF",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("HF",),
            "format_string": "s",
            "description": "Method:",
            "help_text": "The SCF method. HF is the level Bross et al. recommend.",
        },
        "basis": {
            "default": "cc-pVTZ-PP",
            "kind": "string",
            "default_units": "",
            "enumeration": ("cc-pVTZ-PP", "aug-cc-pVTZ-PP", "cc-pVQZ-PP"),
            "format_string": "s",
            "description": "Basis set:",
            "help_text": (
                "The basis set, from the Basis Set Exchange, for the elements it "
                "covers; its pseudopotentials must have published spin-orbit data."
            ),
        },
        "other basis": {
            "default": "cc-pVTZ",
            "kind": "string",
            "default_units": "",
            "enumeration": ("cc-pVTZ", "aug-cc-pVTZ", "cc-pVQZ"),
            "format_string": "s",
            "description": "Basis set for other elements:",
            "help_text": "The basis set for the elements the first one does not cover.",
        },
        "uncontract": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Uncontract the basis sets:",
            "help_text": (
                "Use the basis sets fully uncontracted, as Bross and Peterson "
                "recommend for the spin-orbit term."
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("SpinOrbitParameters.__init__")

        super().__init__(
            defaults={**SpinOrbitParameters.parameters, **defaults}, data=data
        )
