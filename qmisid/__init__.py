"""Columnar (uproot + awkward + numpy + iminuit) electron charge-flip framework.

Steps:
    1. ``truth``       - truth-matched rates in MC        (ChargeFlipTruthv5New.C)
    2. ``likelihood``  - likelihood rates from Z->ee       (Likelihoodv5New.C)
    3. ``closure``     - closure test / rate application  (ClosureTestv5New.C)
    4. ``systematics`` - systematic uncertainties on the rates (new)
"""

__version__ = "0.1.0"
