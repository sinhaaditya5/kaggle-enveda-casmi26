"""Library-matching baseline with a class-1-like validation.

Validation: sample structures seen in >= 2 libraries; their spectra from one library are the queries and
are removed from the reference library. Candidates = train structures within a mass window of the query
neutral mass; score = mean over query spectra of the best cosine to any library spectrum of the candidate.

Usage: python src/baseline.py [n_val] [ppm]
"""
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds
import scipy.sparse as sp

from common import ADDUCT_DELTA, RAW, mrr_at_25, neutral_mass
from structs import load

BIN = 0.01
MAX_PEAKS = 100
PPM_CLEAN = 10.0


def to_sparse(mzs, ints, spread):
    """Rows = spectra. spread=True also puts each query peak in the neighbouring bins (~±BIN tolerance)."""
    rows, cols, vals = [], [], []
    for i, (m, x) in enumerate(zip(mzs, ints)):
        m, x = np.asarray(m), np.asarray(x)
        if len(m) > MAX_PEAKS:
            top = np.argpartition(-x, MAX_PEAKS)[:MAX_PEAKS]
            m, x = m[top], x[top]
        b = np.round(m / BIN).astype(np.int64)
        w = np.sqrt(x)
        for d in ((-1, 0, 1) if spread else (0,)):
            rows.append(np.full(len(b), i)); cols.append(b + d); vals.append(w)
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(mzs), int(3e6)))
    M.sum_duplicates()
    return M


def unit_norm(M, ref=None):
    n = np.sqrt(np.asarray((ref if ref is not None else M).multiply(ref if ref is not None else M).sum(1))).ravel()
    n[n == 0] = 1
    return sp.diags(1 / n) @ M


def main(n_val=500, ppm=5.0, seed=0):
    t0 = time.time()
    S = load()
    mass = S.mass.values
    key2i = {k: i for i, k in enumerate(S.inchikey14)}

    # --- choose validation structures: >=2 libs, mass in test range
    multi = S[(S.libs.str.contains(",")) & S.mass.between(157, 1159)]
    val = multi.sample(n_val, random_state=seed)
    val_keys = set(val.inchikey14)

    # --- load spectra of val structures (queries) first, to pick query lib
    dset = ds.dataset(RAW / "train.parquet")
    cols = ["inchikey14", "ingest_lib", "adduct", "precursor_mz", "precursor_error_ppm", "ms2_mzs",
            "ms2_normalized_intensities"]
    flt = pc.field("adduct").isin(list(ADDUCT_DELTA)) & (pc.abs(pc.field("precursor_error_ppm")) < PPM_CLEAN)
    q_all = dset.to_table(columns=cols, filter=flt & pc.field("inchikey14").isin(list(val_keys))).to_pandas()
    rng = np.random.default_rng(seed)
    qlib = {k: rng.choice(sorted(g.ingest_lib.unique())) for k, g in q_all.groupby("inchikey14")
            if g.ingest_lib.nunique() >= 2}
    q_all = q_all[q_all.inchikey14.isin(qlib)]
    q = q_all[q_all.ingest_lib == q_all.inchikey14.map(qlib)].reset_index(drop=True)
    print(f"val molecules {q.inchikey14.nunique()}, query spectra {len(q)}  [{time.time()-t0:.0f}s]")

    # --- candidates per molecule from mass window
    q["M"] = [neutral_mass(m, a) for m, a in zip(q.precursor_mz, q.adduct)]
    cand = {}
    for k, g in q.groupby("inchikey14"):
        m = g.M.median()
        lo, hi = np.searchsorted(mass, [m * (1 - ppm * 1e-6), m * (1 + ppm * 1e-6)])
        cand[k] = np.arange(lo, hi)
    sizes = np.array([len(v) for v in cand.values()])
    print(f"candidates/molecule: median {np.median(sizes):.0f}, mean {sizes.mean():.0f}, max {sizes.max()}; "
          f"truth in window {np.mean([key2i[k] in set(v) for k, v in cand.items()]):.1%}")
    cand_keys = set(S.inchikey14.values[np.unique(np.concatenate(list(cand.values())))])

    # --- reference library for candidate structures (query lib spectra of that structure removed)
    ref = dset.to_table(columns=cols, filter=flt & pc.field("inchikey14").isin(list(cand_keys))).to_pandas()
    print(f"reference spectra {len(ref):,}  [{time.time()-t0:.0f}s]")
    ref = ref.merge(pd.Series(qlib, name="qlib").rename_axis("inchikey14").reset_index(), how="left")
    ref = ref[~((ref.ingest_lib == ref.qlib))].reset_index(drop=True)  # drop query lib for val structures
    Rm = unit_norm(to_sparse(ref.ms2_mzs, ref.ms2_normalized_intensities, spread=False))
    ref_groups = ref.groupby("inchikey14").indices

    Qm_raw = to_sparse(q.ms2_mzs, q.ms2_normalized_intensities, spread=True)
    Qn = unit_norm(Qm_raw, to_sparse(q.ms2_mzs, q.ms2_normalized_intensities, spread=False))
    sims = (Qn @ Rm.T).tocsr()  # query spectra x reference spectra
    print(f"similarity matrix done  [{time.time()-t0:.0f}s]")

    ranked, truth = [], []
    for k, g in q.groupby("inchikey14"):
        qi = g.index.values
        sub = sims[qi]
        scored = []
        for ci in cand[k]:
            ck = S.inchikey14.values[ci]
            r = ref_groups.get(ck)
            if r is None:
                continue
            best = sub[:, r].max(axis=1).toarray().ravel()  # best match per query spectrum
            scored.append((best.mean(), S.smiles.values[ci], ck))
        scored.sort(reverse=True)
        ranked.append([s[2] for s in scored]); truth.append(k)
    mrr = mrr_at_25(ranked, truth)
    hit = np.mean([k in r[:25] for r, k in zip(ranked, truth)])
    top1 = np.mean([bool(r) and r[0] == k for r, k in zip(ranked, truth)])
    print(f"MRR@25 {mrr:.4f}  top1 {top1:.3f}  recall@25 {hit:.3f}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 500, float(sys.argv[2]) if len(sys.argv) > 2 else 5.0)
