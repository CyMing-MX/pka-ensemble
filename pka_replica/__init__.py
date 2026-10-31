"""pka_replica: a local replica of Rowan's macroscopic/microscopic pKa prediction workflows.

Architecture (mirrors labs.rowansci.com macropka + pka workflows):

    SMILES
      -> microstate enumeration        (SMARTS protonation rules, alternating acid/base rounds)
      -> per-microstate free energies  (ML site-pKa backend + optional electrostatic coupling)
      -> Boltzmann partition function  (Uni-pKa / Starling thermodynamics, eqs. 1-4)
      -> macro-pKa, pI, logD(pH), microstate populations

References:
  - Luo et al., "Bridging Machine Learning and Thermodynamics for Accurate pKa Prediction" (JACS Au, Uni-pKa)
  - Wagen et al., "Physics-Informed Machine Learning Enables Rapid Macroscopic pKa Prediction" (Starling)
  - docs.rowansci.com/workflows/pka, docs.rowansci.com/science/workflows/macropka
"""

__version__ = "0.1.0"

from .workflow import run_macropka, run_micropka  # noqa: F401
