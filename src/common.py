"""Shared helpers: adduct mass table, neutral-mass maths, MRR@25 with the competition's matching rule."""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"

PROTON = 1.007276
H2O = 18.010565
# ion m/z = neutral mass + ADDUCT_DELTA[adduct]
ADDUCT_DELTA = {
    "[M+H]+": PROTON,
    "[M+NH4]+": 18.033823,
    "[M+Na]+": 22.989218,
    "[M+K]+": 38.963158,
    "[M-H2O+H]+": PROTON - H2O,
    "[M-2H2O+H]+": PROTON - 2 * H2O,
    "[M-H]-": -PROTON,
    "[M-H2O-H]-": -PROTON - H2O,
    "[M+CH2O2-H]-": 46.005479 - PROTON,
    "[M+Cl]-": 34.969402 + 0.000549,
}


def neutral_mass(precursor_mz, adduct):
    return precursor_mz - ADDUCT_DELTA[adduct]


def inchikey14(smiles):
    """First InChIKey block after RDKit tautomer canonicalisation (the competition's matching key)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import inchi
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.TautomerEnumerator().Canonicalize(mol)
        return inchi.MolToInchiKey(mol)[:14]
    except Exception:
        return None


def mrr_at_25(ranked, truth):
    """ranked: list of lists of keys (best first); truth: list of keys."""
    s = 0.0
    for r, t in zip(ranked, truth):
        for i, k in enumerate(r[:25], 1):
            if k == t:
                s += 1.0 / i
                break
    return s / len(truth)
