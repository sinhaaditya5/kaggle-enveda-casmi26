"""Metadata EDA for Enveda CASMI 2026. Reads only scalar columns of train (peak arrays are skipped).

Usage: python src/eda.py  ->  writes reports/eda.md
"""
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
OUT = Path(__file__).resolve().parent.parent / "reports" / "eda.md"

SCALAR = ["base_peak_intensity", "adduct", "ionization_mode",
          "instrument_type", "precursor_mz", "collision_energy_ev", "normalized_smiles",
          "inchikey14", "molecular_formula", "ingest_lib", "precursor_error_ppm", "num_peaks"]


def md(df, **kw):
    return df.to_markdown(**kw)


def main():
    tr = pq.read_table(RAW / "train.parquet", columns=SCALAR).to_pandas()
    te = pq.read_table(RAW / "test.parquet", columns=["molecule_id", "spectrum_id"] + [
        c for c in SCALAR if c in pq.read_schema(RAW / "test.parquet").names]).to_pandas()
    tr["n_ce"] = tr.collision_energy_ev.map(lambda x: 0 if x is None else len(x))
    te["n_ce"] = te.collision_energy_ev.map(len)
    test_adducts = set(te.adduct.unique())
    out = ["# EDA\n", f"Train: {len(tr):,} spectra, {tr.inchikey14.nunique():,} structures. "
           f"Test sample: {len(te):,} spectra, {te.molecule_id.nunique():,} molecules.\n"]

    lib = tr.groupby("ingest_lib").agg(spectra=("adduct", "size"), structures=("inchikey14", "nunique"),
                                       med_peaks=("num_peaks", "median"),
                                       pos_frac=("ionization_mode", lambda s: (s == "positive").mean()),
                                       med_ppm=("precursor_error_ppm", lambda s: s.abs().median()))
    out += ["## Libraries\n", md(lib.sort_values("spectra", ascending=False), floatfmt=".2f"), "\n"]

    ad = tr.adduct.value_counts().head(25).to_frame("train")
    ad["in_test_set"] = ad.index.isin(test_adducts)
    out += ["## Top train adducts\n", md(ad), "\n",
            f"Share of train spectra with a test-set adduct: {tr.adduct.isin(test_adducts).mean():.1%}\n"]

    out += ["## Test adducts / mode\n", md(te.adduct.value_counts().to_frame("spectra")), "\n",
            md(te.ionization_mode.value_counts().to_frame("spectra")), "\n"]

    out += ["## Spectra per structure (train)\n", md(tr.groupby("inchikey14").size().describe().to_frame("n_spectra")), "\n",
            "## Spectra per molecule (test sample)\n", md(te.groupby("molecule_id").size().describe().to_frame("n_spectra")), "\n"]

    out += ["## Collision-energy list length\n", md(pd.DataFrame({"train": tr.n_ce.value_counts(),
                                                               "test": te.n_ce.value_counts()}).fillna(0).astype(int)), "\n"]

    out += ["## Precursor error (ppm), train\n", md(tr.precursor_error_ppm.abs().describe(
        percentiles=[.5, .9, .99]).to_frame("abs_ppm"), floatfmt=".2f"), "\n"]

    out += ["## Precursor m/z, train vs test\n", md(pd.DataFrame({"train": tr.precursor_mz.describe(),
                                                                "test": te.precursor_mz.describe()}), floatfmt=".1f"), "\n"]

    inst = tr.loc[tr.ingest_lib.isin(["enveda-180", "enveda-np-examples"])]
    out += [f"## Instrument-matched rows (enveda-*)\n\n{len(inst):,} spectra, {inst.inchikey14.nunique():,} structures.\n"]

    OUT.write_text("\n".join(out))
    print(OUT)


if __name__ == "__main__":
    main()
