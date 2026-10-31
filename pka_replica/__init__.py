"""pka_replica: aqueous microscopic & macroscopic pKa prediction.

Pipeline:

    SMILES
      -> microstate enumeration        (SMARTS protonation rules, alternating acid/base rounds)
      -> per-microstate free energies  (fragment-constant backend, or GFN2-xTB/ALPB physics
                                       with optional AIMNet2 pre-optimization)
      -> Boltzmann partition function  (Uni-pKa / Starling thermodynamics, eqs. 1-4)
      -> macro-pKa, pI, logD(pH), microstate populations

References:
  - Luo et al., "Bridging Machine Learning and Thermodynamics for Accurate
    pKa Prediction" (JACS Au, Uni-pKa)
  - Wagen et al., "Physics-Informed Machine Learning Enables Rapid
    Macroscopic pKa Prediction" (Starling, 2025)
  - Bannwarth et al., GFN2-xTB and ALPB (xtb)
  - Sipos-Szabó et al., pKahub benchmark (JCIM 2026)
"""

__version__ = "1.0.0"

from .workflow import run_macropka, run_micropka  # noqa: F401
