# DiRT Candidate Index Gene Finder

This repository provides a reference implementation for selecting **candidate index genes** per target gene using the DiRT concept (DEG-by-index Ratio Transformation). For each target gene, the script ranks all other genes by the **normalized dispersion** (NDIV = std/mean) of control-only ratios and returns the top `N` candidates that best preserve baseline co‑expression with the target.

> **TL;DR**: Given a CPM matrix of 10,000 genes × 18 samples (C1–C9 controls, T1–T9 treated), the script computes, for each target gene, the 10 candidates with the lowest dispersion of control ratios and outputs their across‑all‑samples ratios together with the NDIV score.

---

## Contents

- `script.py` (your working code block)
- `README.md` (this file)

---

## Input

- **File**: `Test.csv`
- **Required columns**:
  - `Geneid` (unique gene identifier per row)
  - Control samples: `C1 ... C9`
  - Treated samples: `T1 ... T9`
- **Values**: CPM (floating‑point numeric). Missing values should be imputed or removed prior to running.

**Example header**:

```
Geneid,C1,C2,C3,C4,C5,C6,C7,C8,C9,T1,T2,T3,T4,T5,T6,T7,T8,T9
```

> If your column names differ or you have a different number of samples, adjust the `control_start`, `control_end`, `all_start`, and `all_end` parameters in `tencand(...)` accordingly.

---

## Output

- **File**: `DiRT_test.csv` (or `DiRT_test1.csv`, `DiRT_test2.csv` for split runs)
- **Rows**: Concatenated top‑`N` candidate rows for each processed target index.
- **Columns**:
  - `ID` — string in the form `targetGene/candidateGene`
  - `ndiv` — normalized dispersion (std/mean) of control ratios for that pair
  - Ratio columns for **all** samples (`C1..C9`, `T1..T9`) = `(target / candidate)`

**Preview** (first few columns):

```
ID,ndiv,C1,C2,...,C9,T1,...,T9
ATG1/ACTB,0.0421,1.03,0.98,...
...
```

---

## How it works

For a given target gene at row index `idx`:

1. Take the target’s control vector `target[C1:C9]`.
2. For each candidate gene `x`, compute control ratios `target[C1:C9] / x[C1:C9]`.
3. Compute **NDIV** = `std(ratios, ddof=1) / mean(ratios)`.
4. Rank all candidates by NDIV ascending, **excluding the self‑pair**.
5. Return the top `N` candidates, with their **full** ratios across `C1..C9,T1..T9` and the NDIV score.

A small epsilon (`EPS = 1e-12`) is added to denominators to avoid divide‑by‑zero.

---

## Usage

1. Install dependencies (Python ≥ 3.9 recommended):

```bash
pip install pandas numpy
```

2. Place `Test.csv` in the working directory.

3. Run the code block (for a quick test it processes the first 10 targets and writes `DiRT_test.csv`).

**Modifying ranges**

- To process the full matrix (10,000 target indices), replace the test loop with the “full run” loops in the comments, or parallelize (see below).

**Function parameters**

- `tencand(idx, control_start="C1", control_end="C9", all_start="C1", all_end="T9", top_n=10)`
  - `idx`: integer row index of the target gene.
  - `control_start`, `control_end`: inclusive labels defining the control block.
  - `all_start`, `all_end`: inclusive labels defining the full (control+treated) block.
  - `top_n`: number of candidate index genes to return (default 10).

---

## Performance & Scaling

This implementation calls `DataFrame.apply` over **all genes** for each target gene, making the runtime roughly **O(G²)** for `G` genes. For 10,000 genes, that’s heavy.

**Practical tips**

- **Split the workload** (already shown in the script):
  - Half 1: targets `[0, 5000)` → `DiRT_test1.csv`
  - Half 2: targets `[5000, 10000)` → `DiRT_test2.csv`
  - Concatenate results.
- **Persist intermediate chunks** every few hundred targets to avoid data loss.
- **Use Feather/Parquet** for fast IO if you iterate multiple times.

**Optional parallelization (sketch)**

```python
from multiprocessing import Pool, cpu_count

def run_block(start, stop):
    out = pd.concat([tencand(i) for i in range(start, stop)], axis=0, ignore_index=True)
    out.to_csv(f"DiRT_{start}_{stop}.csv", index=False)

if __name__ == "__main__":
    G = len(df)
    blocks = [(0, 5000), (5000, 10000)]  # adjust
    with Pool(processes=min(len(blocks), cpu_count())) as p:
        p.starmap(run_block, blocks)
```

**Roadmap for speed‑ups** (future work)

- Vectorize NDIV using matrix operations (precompute log‑ratios and use variance identities).
- Use **Numba** or **CuPy** for JIT/GPU acceleration.
- Batch computations to reuse temporary arrays and reduce Python overhead.

---

## Reproducibility

- Deterministic given identical `Test.csv` and environment.
- Set pandas options to ensure consistent float formatting if needed.

---

## Troubleshooting

- **ValueError: cannot do slice indexing with these indexers**  
  Column labels must exist and be **ordered** so that `control_start:control_end` and `all_start:all_end` form **contiguous** blocks.
- **NaNs or infs in output**  
  Ensure all CPM values are numeric and non‑negative; impute or filter zeros if ratios become unstable (EPS reduces but cannot fully prevent issues with all‑zero rows).
- **Duplicate `Geneid`**  
  Recommended that `Geneid` is unique per row to avoid ambiguity in `ID` strings.

---

## Method notes (NDIV)

- **NDIV = std/mean** of control ratios quantifies **relative dispersion**. Lower NDIV means the candidate gene’s control expression tracks the target gene more closely (stable baseline relationship), which is desirable for forming a ratio‑based index.
- Self‑pairs are removed by testing the `ID` suffix `"/{target_geneid}"`.

---

## License

MIT (or your preferred license).

---

## Citation

If you use this script in a publication, please cite your DiRT paper or preprint and acknowledge this repository.
