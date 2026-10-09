# Enveda CASMI 2026 · Molecule ID from Mass Spectra

My work on the [Enveda CASMI 2026](https://www.kaggle.com/competitions/enveda-CASMI26-molecule-id-mass-spectra) Kaggle featured code competition:
from tandem mass spectra (LC-MS/MS) of a natural-product-like molecule, predict its 2D structure as SMILES.

| | |
|---|---|
| **Task** | Per molecule, a ranked list of up to 25 candidate SMILES (best first, `;`-joined) |
| **Metric** | MRR@25: 1/rank of the first correct structure, 0 if absent. Match = same InChIKey14 after RDKit tautomer canonicalisation (stereo ignored) |
| **Train** | `train.parquet`: ~2.5 M spectra, ~275 k structures, 9 source libraries (`ingest_lib`), ~3 GB |
| **Test** | ~1,500 spectra of ~400 molecules (1-16 spectra each, median 3), all Bruker timsTOF, 157-1,159 Da. Hidden set replaces the local sample at re-run |
| **Novelty classes** | 1: in public spectral libraries · 2: known structure in PubChem/COCONUT, no public spectra · 3: novel, needs de novo generation |
| **Submission** | Kaggle notebook, CPU or GPU ≤ 9 h, no internet, `submission.csv`; public external data and pretrained models allowed |
| **Timeline** | 14 Sep 2026 → entry/merger deadline 7 Dec → **final 14 Dec 2026** (23:59 UTC) |
| **Prize** | $50,000 · medals and points |

## The core problem

Library matching solves class 1 only. Classes 2 and 3 need a model that maps a spectrum to a structure,
so the pipeline is likely three parts:

1. **Aggregate** all spectra of a molecule (different adducts and collision energies) into one representation.
2. **Retrieve**: embed spectrum and candidate structures (PubChem / COCONUT / train structures) in a shared space, or predict a fingerprint and rank candidates by formula + fingerprint similarity.
3. **Generate** (class 3): de novo candidates conditioned on formula and predicted fingerprint, then re-rank.

Formula is strongly constrained by the precursor m/z and adduct, which prunes any candidate list hard.

## Key findings from EDA

Full report: [reports/eda.md](reports/eda.md) (regenerate with `python src/eda.py`).

- **Train is 2.54 M spectra / 275,810 structures** (median 6 spectra per structure). It has **no `molecule_id`/`spectrum_id`**, so group by `inchikey14`.
- **Instrument-matched data is only `enveda-180` + `enveda-np-examples`** (1.15 M spectra, 183 k structures); the rest is Orbitrap/QTOF from other labs.
- **Test adducts are 7 of 26+ train adducts** (84.7 % of train spectra share a test adduct). The local test sample is 79 % `[M+H]+`, 16 % `[M-H]-`.
- **Collision energies**: test is single-CE (75 %) or 3-CE merged (25 %); train also has 337 k spectra with no CE at all.
- **Precursor error is dirty**: median 1.5 ppm but the 99th percentile is ~2,900 ppm, so filter train on `precursor_error_ppm`.
- **The local `test.parquet` is a sample of train** (1,213 spectra, 400 molecules), so it cannot be used for validation. The hidden test has a different distribution (precursor m/z 245-460 here vs up to 1,159 Da hidden).

## Plan

| Stage | What | Status |
|---|---|---|
| 0 | Accept rules, download data, EDA (library mix, adducts, CE, overlap of train structures with test class 1) | pending |
| 1 | Baseline: cosine library match against train spectra, per-molecule aggregation, MRR@25 validation split by structure | pending |
| 2 | Fingerprint / embedding model + candidate retrieval from public structure databases (offline, as a Kaggle dataset) | pending |
| 3 | Formula-constrained de novo generation for class 3 | pending |
| 4 | Ensemble and re-rank, two final picks | pending |

## Layout

`src/` code · `reports/` EDA and experiment notes · `notebooks/` Kaggle kernels · `docs/` figures and the draw.io pipeline.
Data is never committed (Kaggle rules); `data/` is git-ignored.
