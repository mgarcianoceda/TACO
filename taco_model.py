"""TACO model — self-contained scorer for TCR-pMHC activation prediction.

Load a fitted model from its JSON spec and score TCR-peptide pairs:

    from taco_model import TACOModel

    m = TACOModel.load("taco_model.json")
    df["score"] = m.score(df)
    df["prob"]  = m.predict_proba(df)
    panel = m.score(df_tcr, peptides=["GILG...", "NLVP..."])
    df["score_beta"] = m.score_beta(df)

Requires only numpy and pandas (no sklearn, no torch).
"""

import json
import numpy as np
import pandas as pd

class TACOModel:
    """A fitted TACO model, loaded from a file, scoring dataframes.

    SELF-CONTAINED ON PURPOSE.  This class rebuilds the contact geometry from the amino-acid sequences and the constants in its own file.
    Copy the class and the file together and inference works anywhere numpy and pandas do.

        m = TACOModel.load("path/to/taco_model.json")

        # Per-row scoring — df has the six CDR columns and `peptide`
        df["score"] = m.score(df)
        df["score"] = m.score(df, parallel=True)

        # Panel scoring — every TCR scored against every peptide, no cross-product DataFrame
        panel = m.score(df_tcr, peptides=["GILG...", "NLVP..."])  # returns DataFrame

        # Beta-only scoring — alpha loops excluded
        df["score_beta"] = m.score_beta(df)
        panel_beta = m.score_beta(df_tcr, peptides=["GILG...", "NLVP..."])

    The score is a transferable log-odds with a single global intercept, so it is directly
    comparable across antigens — no normalization needed."""


    COLS = ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b", "peptide"]
    _ALPHA_LOOPS = {"cdr1a", "cdr2a", "cdr3a"}

    def __init__(self, spec):
        self.spec = spec
        self.aa = {c: i for i, c in enumerate(spec["alphabet"])}

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            spec = json.load(fh)
        if spec.get("format") != "taco-model-1":
            raise ValueError(f"{path}: not an TACO model file")
        return cls(spec)

    # ---- geometry ------------------------------------------------------------------------
    @staticmethod
    def _w(w9, Lp):
        """Nine control points -> one weight per peptide position."""
        w9 = np.asarray(w9, np.float32)
        if Lp == 9:
            return w9
        x9 = (np.arange(1, 10) - 0.5) / 9
        return np.interp((np.arange(1, Lp + 1) - 0.5) / Lp, x9, w9).astype(np.float32)

    def _idx(self, seqs, Lmax):
        aa = self.aa
        n = len(seqs)
        I = np.full((n, Lmax), -1, np.int16)
        for r in range(n):
            s = str(seqs[r])
            end = min(len(s), Lmax)
            for j in range(end):
                I[r, j] = aa.get(s[j], -1)
        return I

    def _loop(self, I, L, geo, Lp):
        """Contact geometry of one loop against a peptide of length Lp: (n_tcr, Lp, 20)."""
        n, Lmax = I.shape
        T = np.zeros((n, Lp, 20), np.float32)
        i = np.arange(1, Lmax + 1, dtype=np.float32)[None, :]
        j = np.arange(1, Lp + 1, dtype=np.float32)[None, None, :]
        Lb = L.astype(np.float32)[:, None]
        ctr = (Lb + 1) / 2
        u = np.exp(-self.spec["kappa"] * (i - ctr) ** 2 / Lb) * (i <= Lb)
        s = (geo["centre"] + geo["slope"] * (i - ctr)) * (Lp / 9.0)
        C = (u[:, :, None] * np.exp(-((s[:, :, None] - j) ** 2) / (2 * geo["sigma"] ** 2))).astype(
            np.float32
        )
        for l in range(Lmax):
            col = I[:, l]
            for a in range(20):
                m = col == a
                if m.any():
                    T[m, :, a] += C[m, l, :]
        return T

    # ---- TCR dedup + encoding (shared across all paths) ---------------------------------
    def _dedup_and_encode(self, df_tcr):
        """Integer-hash dedup, encode sequences, compute geometry-ready arrays.

        Returns (codes, IDX, LEN, Lmax) where codes maps each input row to a unique-TCR index."""
        loops = [c for c in self.COLS if c != "peptide"]
        combined = np.zeros(len(df_tcr), dtype=np.int64)
        for c in loops:
            fc, _ = pd.factorize(df_tcr[c].values)
            combined = combined * 100003 + fc
        codes, _ = pd.factorize(combined)
        uniq_idx = np.unique(codes, return_index=True)[1]

        SEQ = {c: df_tcr[c].values[uniq_idx].astype(str) for c in loops}
        LEN = {c: np.array([len(s) for s in SEQ[c]]) for c in loops}
        Lmax = int(max(v.max() for v in LEN.values()))
        IDX = {c: self._idx(SEQ[c], Lmax) for c in loops}
        return codes, IDX, LEN, Lmax

    def _compute_geometry(self, IDX, LEN, Lmax, peptide_lengths):
        """Compute surface geometry for each unique peptide length."""
        loops = [c for c in self.COLS if c != "peptide"]
        all_S = {}
        for Lp in sorted(set(peptide_lengths)):
            all_S[Lp] = {
                su["name"]: sum(
                    self._loop(IDX[c], LEN[c], self.spec["loops"][c], Lp)
                    for c in su["loops"]
                )
                for su in self.spec["surfaces"]
            }
        return all_S

    # ---- setup for DataFrame mode -------------------------------------------------------
    @staticmethod
    def _build_pep_index(pep):
        """Pre-compute {length: {peptide_str: row_indices}} via numpy grouping."""
        unique_peps, inverse = np.unique(pep, return_inverse=True)
        order = np.argsort(inverse, kind='mergesort')
        counts = np.bincount(inverse)
        splits = np.concatenate([[0], np.cumsum(counts)])
        pep_idx = {}
        for k, p in enumerate(unique_peps):
            Lp = len(p)
            if Lp not in pep_idx:
                pep_idx[Lp] = {}
            pep_idx[Lp][p] = order[splits[k]:splits[k + 1]]
        return pep_idx

    def _score_setup(self, df):
        """Dedup, encode, geometry, peptide index — everything reusable across peptides."""
        codes, IDX, LEN, Lmax = self._dedup_and_encode(df)
        pep = df["peptide"].values.astype(str)
        pep_idx = self._build_pep_index(pep)
        all_S = self._compute_geometry(IDX, LEN, Lmax, [Lp for Lp in pep_idx])
        return codes, pep_idx, all_S

    # ---- the score -----------------------------------------------------------------------
    def score(self, df, peptides=None, parallel=False, n_workers=None):
        """Score TCR-peptide pairs.

        Basic:     df["score"] = m.score(df)
        Panel:     panel = m.score(df_tcr, peptides=["AA","BB"])   # DataFrame (n_tcr x n_pep)

        Scores include a global bias and are directly comparable across antigens."""
        loops = [c for c in self.COLS if c != "peptide"]

        # --- Panel mode (peptides given explicitly) ---
        if peptides is not None:
            peptides = list(peptides)
            missing = [c for c in loops if c not in df.columns]
            if missing:
                raise ValueError(f"dataframe is missing {missing}")
            panel = self._score_panel(df, peptides, parallel, n_workers)
            return pd.DataFrame(panel, index=df.index, columns=peptides)

        # --- DataFrame mode (df has peptide column) ---
        missing = [c for c in self.COLS if c not in df.columns]
        if missing:
            raise ValueError(
                f"dataframe is missing {missing}; TACOModel.score needs "
                f"{self.COLS} (upper-case one-letter sequences, no gaps)"
            )
        if parallel:
            return self._score_parallel(df, n_workers)
        return self._score_single(df)

    def predict_proba(self, df, peptides=None, parallel=False, n_workers=None):
        """Predicted activation probability (sigmoid of score).

        Same interface as score(). Returns values in [0, 1]."""
        s = self.score(df, peptides=peptides, parallel=parallel, n_workers=n_workers)
        return 1.0 / (1.0 + np.exp(-s))

    def _compute_geometry_beta(self, IDX, LEN, Lmax, peptide_lengths):
        """Geometry excluding alpha loops — only beta chain contributions."""
        all_S = {}
        for Lp in sorted(set(peptide_lengths)):
            all_S[Lp] = {}
            for su in self.spec["surfaces"]:
                beta_loops = [c for c in su["loops"] if c not in self._ALPHA_LOOPS]
                if beta_loops:
                    all_S[Lp][su["name"]] = sum(
                        self._loop(IDX[c], LEN[c], self.spec["loops"][c], Lp)
                        for c in beta_loops
                    )
                else:
                    n_uniq = next(iter(IDX.values())).shape[0]
                    all_S[Lp][su["name"]] = np.zeros((n_uniq, Lp, 20), np.float32)
        return all_S

    def score_beta(self, df, peptides=None, parallel=False, n_workers=None):
        """Score using only TCR-beta chains; alpha loop contributions are excluded entirely.

        Basic:     df["score_beta"] = m.score_beta(df)
        Panel:     panel = m.score_beta(df_tcr, peptides=["AA","BB"])

        Alpha loops are dropped from the geometry sum rather than replaced with dummies.
        The score retains the global bias and cross-antigen comparability — no per-peptide
        normalization is applied."""
        if peptides is not None:
            return self._score_beta_panel(df, peptides)
        return self._score_beta_single(df)

    def _score_beta_single(self, df):
        codes, IDX, LEN, Lmax = self._dedup_and_encode(df)
        pep = df["peptide"].values.astype(str)
        pep_idx = self._build_pep_index(pep)
        all_S = self._compute_geometry_beta(IDX, LEN, Lmax, [Lp for Lp in pep_idx])
        return self._score_from_setup(codes, pep_idx, all_S, len(df))

    def _score_beta_panel(self, df_tcr, peptides):
        codes, IDX, LEN, Lmax = self._dedup_and_encode(df_tcr)
        pep_lengths = [len(p) for p in peptides]
        all_S = self._compute_geometry_beta(IDX, LEN, Lmax, pep_lengths)
        aa = self.aa
        n_unique = int(codes.max()) + 1 if len(codes) > 0 else 0
        unique_out = np.zeros((n_unique, len(peptides)))
        for pi, p in enumerate(peptides):
            Lp = len(p)
            b = np.array([aa.get(ch, -1) for ch in p])
            if (b < 0).any():
                raise ValueError("peptide has a residue outside the alphabet")
            z = np.zeros(n_unique)
            for su in self.spec["surfaces"]:
                M = np.asarray(su["M"], np.float32)
                w = self._w(su["w9"], Lp)
                G = all_S[Lp][su["name"]]
                for j in range(Lp):
                    z += (G[:, j, :] @ M[:, b[j]]) * w[j]
            unique_out[:, pi] = z
        return pd.DataFrame(unique_out[codes] + self.spec.get("bias", 0.0),
                            index=df_tcr.index, columns=peptides)

    # ---- DataFrame scoring (original path) -----------------------------------------------
    def _score_single(self, df):
        """Single-threaded scoring."""
        codes, pep_idx, all_S = self._score_setup(df)
        return self._score_from_setup(codes, pep_idx, all_S, len(df))

    def _score_from_setup(self, codes, pep_idx, all_S, n_rows):
        """Score peptides from precomputed setup. Used by single and parallel paths."""
        aa = self.aa
        out = np.zeros(n_rows)
        for Lp in sorted(pep_idx):
            S = all_S[Lp]
            for p in sorted(pep_idx[Lp]):
                rows = pep_idx[Lp][p]
                t = codes[rows]
                b = np.array([aa.get(ch, -1) for ch in p])
                if (b < 0).any():
                    raise ValueError(f"peptide has a residue outside the alphabet")
                z = np.zeros(len(rows))
                for su in self.spec["surfaces"]:
                    M = np.asarray(su["M"], np.float32)
                    w = self._w(su["w9"], Lp)
                    G = S[su["name"]][t]
                    for j in range(Lp):
                        z += (G[:, j, :] @ M[:, b[j]]) * w[j]
                out[rows] = z
        return out + self.spec.get("bias", 0.0)

    def _score_parallel(self, df, n_workers=None):
        """Parallel scoring with automatic strategy selection."""
        import multiprocessing
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        n_workers = min(n_workers, max(1, len(df) // 500))
        if n_workers <= 1:
            return self._score_single(df)

        pep = df["peptide"].values.astype(str)
        n_unique_pep = len(np.unique(pep))

        if n_unique_pep >= n_workers:
            return self._score_parallel_by_peptide(df, n_workers)
        return self._score_parallel_by_rows(df, n_workers)

    def _score_parallel_by_peptide(self, df, n_workers):
        """Geometry once in parent, distribute peptide scoring to workers."""
        import multiprocessing
        codes, pep_idx, all_S = self._score_setup(df)

        jobs = [(Lp, p) for Lp in sorted(pep_idx) for p in sorted(pep_idx[Lp])]

        global _taco_par_codes, _taco_par_pep_idx, _taco_par_all_S
        global _taco_par_surfaces, _taco_par_aa
        _taco_par_codes = codes
        _taco_par_pep_idx = pep_idx
        _taco_par_all_S = all_S
        _taco_par_surfaces = self.spec["surfaces"]
        _taco_par_aa = self.aa

        batches = [[] for _ in range(n_workers)]
        for i, job in enumerate(jobs):
            batches[i % n_workers].append(job)

        with multiprocessing.Pool(n_workers) as pool:
            results = pool.map(_taco_score_peptide_batch, batches)

        _taco_par_codes = None
        _taco_par_pep_idx = None
        _taco_par_all_S = None
        _taco_par_surfaces = None
        _taco_par_aa = None

        out = np.zeros(len(df))
        for batch_result in results:
            for rows, z in batch_result:
                out[rows] = z
        out += self.spec.get("bias", 0.0)
        return out

    def _score_parallel_by_rows(self, df, n_workers):
        """Split rows across workers — each does full setup + geometry on its chunk."""
        import multiprocessing
        global _taco_par_spec, _taco_par_df
        _taco_par_spec = self.spec
        _taco_par_df = df.reset_index(drop=True)

        chunk = len(_taco_par_df) // n_workers
        ranges = [(i * chunk, (i + 1) * chunk if i < n_workers - 1 else len(_taco_par_df))
                  for i in range(n_workers)]

        with multiprocessing.Pool(n_workers) as pool:
            results = pool.map(_taco_score_row_range, ranges)

        _taco_par_spec = None
        _taco_par_df = None
        return np.concatenate(results)

    # ---- Panel scoring (TCRs x peptides, no cross-product DataFrame) ---------------------
    def _score_panel(self, df_tcr, peptides, parallel=False, n_workers=None):
        """Score every TCR against every peptide. Returns (n_tcr, n_pep) numpy array."""
        if parallel:
            return self._score_panel_parallel(df_tcr, peptides, n_workers)
        return self._score_panel_single(df_tcr, peptides)

    def _panel_setup(self, df_tcr, peptides):
        """Dedup, encode, geometry for panel mode (TCR df + explicit peptide list)."""
        codes, IDX, LEN, Lmax = self._dedup_and_encode(df_tcr)
        pep_lengths = [len(p) for p in peptides]
        all_S = self._compute_geometry(IDX, LEN, Lmax, pep_lengths)
        return codes, all_S

    def _score_panel_single(self, df_tcr, peptides):
        """Score full panel single-threaded. Geometry computed once, reused across peptides."""
        codes, all_S = self._panel_setup(df_tcr, peptides)
        aa = self.aa
        n_unique = int(codes.max()) + 1 if len(codes) > 0 else 0
        unique_out = np.zeros((n_unique, len(peptides)))

        for pi, p in enumerate(peptides):
            Lp = len(p)
            b = np.array([aa.get(ch, -1) for ch in p])
            if (b < 0).any():
                raise ValueError(f"peptide has a residue outside the alphabet")
            z = np.zeros(n_unique)
            for su in self.spec["surfaces"]:
                M = np.asarray(su["M"], np.float32)
                w = self._w(su["w9"], Lp)
                G = all_S[Lp][su["name"]]
                for j in range(Lp):
                    z += (G[:, j, :] @ M[:, b[j]]) * w[j]
            unique_out[:, pi] = z

        return unique_out[codes] + self.spec.get("bias", 0.0)

    def _score_panel_parallel(self, df_tcr, peptides, n_workers=None):
        """Score full panel, parallel across peptides."""
        import multiprocessing
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()
        n_workers = min(n_workers, len(peptides))
        if n_workers <= 1:
            return self._score_panel_single(df_tcr, peptides)

        codes, all_S = self._panel_setup(df_tcr, peptides)
        n_unique = int(codes.max()) + 1 if len(codes) > 0 else 0

        global _taco_par_all_S, _taco_par_surfaces, _taco_par_aa, _taco_par_n_unique
        _taco_par_all_S = all_S
        _taco_par_surfaces = self.spec["surfaces"]
        _taco_par_aa = self.aa
        _taco_par_n_unique = n_unique

        batches = [[] for _ in range(n_workers)]
        for i, p in enumerate(peptides):
            batches[i % n_workers].append((i, p))

        with multiprocessing.Pool(n_workers) as pool:
            results = pool.map(_taco_score_panel_batch, batches)

        _taco_par_all_S = None
        _taco_par_surfaces = None
        _taco_par_aa = None
        _taco_par_n_unique = None

        unique_out = np.zeros((n_unique, len(peptides)))
        for batch_result in results:
            for pi, z in batch_result:
                unique_out[:, pi] = z

        return unique_out[codes] + self.spec.get("bias", 0.0)

# ---- Module-level globals and workers for multiprocessing (fork-inherited) ----------------
_taco_par_codes = None
_taco_par_pep_idx = None
_taco_par_all_S = None
_taco_par_surfaces = None
_taco_par_aa = None
_taco_par_spec = None
_taco_par_df = None
_taco_par_n_unique = None

def _taco_score_peptide_batch(batch):
    """Worker: score a batch of peptides from precomputed geometry in parent globals."""
    codes = _taco_par_codes
    pep_idx = _taco_par_pep_idx
    all_S = _taco_par_all_S
    surfaces = _taco_par_surfaces
    aa = _taco_par_aa
    results = []
    for Lp, p in batch:
        rows = pep_idx[Lp][p]
        t = codes[rows]
        b = np.array([aa.get(ch, -1) for ch in p])
        z = np.zeros(len(rows))
        for su in surfaces:
            M = np.asarray(su["M"], np.float32)
            w = TACOModel._w(su["w9"], Lp)
            G = all_S[Lp][su["name"]][t]
            for j in range(Lp):
                z += (G[:, j, :] @ M[:, b[j]]) * w[j]
        results.append((rows, z))
    return results

def _taco_score_row_range(start_end):
    """Worker: full scoring on a row slice from fork-inherited globals."""
    start, end = start_end
    m = TACOModel(_taco_par_spec)
    return m._score_single(_taco_par_df.iloc[start:end].reset_index(drop=True))

def _taco_score_panel_batch(batch):
    """Worker: score a batch of peptides against all unique TCRs (panel mode)."""
    all_S = _taco_par_all_S
    surfaces = _taco_par_surfaces
    aa = _taco_par_aa
    n_unique = _taco_par_n_unique
    results = []
    for pi, p in batch:
        Lp = len(p)
        b = np.array([aa.get(ch, -1) for ch in p])
        z = np.zeros(n_unique)
        for su in surfaces:
            M = np.asarray(su["M"], np.float32)
            w = TACOModel._w(su["w9"], Lp)
            G = all_S[Lp][su["name"]]
            for j in range(Lp):
                z += (G[:, j, :] @ M[:, b[j]]) * w[j]
        results.append((pi, z))
    return results
