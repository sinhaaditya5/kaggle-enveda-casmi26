"""Build the structure table (one row per inchikey14: smiles, exact mass, #libs) once and cache it."""
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from common import RAW, ROOT

CACHE = ROOT / "data" / "cache" / "structs.parquet"


def build():
    RDLogger.DisableLog("rdApp.*")
    df = pq.read_table(RAW / "train.parquet",
                       columns=["inchikey14", "normalized_smiles", "ingest_lib"]).to_pandas()
    g = df.groupby("inchikey14").agg(smiles=("normalized_smiles", "first"), n_spec=("ingest_lib", "size"),
                                     libs=("ingest_lib", lambda s: ",".join(sorted(set(s))))).reset_index()
    mass = []
    for s in g.smiles:
        m = Chem.MolFromSmiles(s)
        mass.append(Descriptors.ExactMolWt(m) if m is not None else np.nan)
    g["mass"] = mass
    g = g.dropna(subset=["mass"]).sort_values("mass").reset_index(drop=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    g.to_parquet(CACHE)
    return g


def load():
    return pd.read_parquet(CACHE) if CACHE.exists() else build()


if __name__ == "__main__":
    g = load()
    print(len(g), g.mass.describe())
