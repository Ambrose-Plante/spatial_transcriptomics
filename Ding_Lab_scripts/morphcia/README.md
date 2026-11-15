# MorphCIA — discovering spatial gene/cell-type logic associations

## Research objective

`morphcia.py` addresses the second question in the accompanying `logic_search.pdf` presentation: **“How does [the] donor cell influence [the] receiver cell?”** (slides 23–39). The working hypothesis is that, among receiver cells located near a chosen donor population, the expression status of a **receiver gene** can be associated with the composition and gene expression of the **surrounding cellular neighborhood**.

The objective is **candidate association discovery**, not proof of causal cellular influence. Slide 37 illustrates an AND-rule candidate involving `fibroblast_WNT5A / C1R` and `tumor / IQGAP2`; such findings require independent biological and computational validation.

## Implemented algorithm in `morphcia.py`

1. **Select sections:** read `sample_organs.txt`, enumerate section datasets, and optionally restrict to `colon` or `liver` (`all` means all available sections). Some colon/liver section pairs are manually listed in the script.
2. **Choose cell classes:** select receiver cells (configured as tumor classes) and a donor population (configured as `fibroblast_WNT5A`). Use a `scipy.spatial.cKDTree` to find the nearest donor for each receiver and retain receiver cells within `dist_thresh`.
3. **Condition on receiver-gene expression:** for each configured receiver gene, split those near-donor receivers into **positive** (expression >0) and **negative** (expression ==0) groups.
4. **Build environment representations:** for each positive or negative receiver cell, inspect cells of *all recorded cell types* within `dist_thresh`. For every cell-type × gene combination, aggregate expression using the configured `max`, `sum` or `mean`. The resulting arrays are effectively **receiver cell × cell type × gene** matrices, conceptually corresponding to the presentation's “donoractive” / “donorinactive” feature arrays (slides 27–30).
5. **Pre-filter candidate features:** retain a gene/cell-type feature if it is expressed above zero in **>30%** of at least one of the two receiver groups.
6. **Record single-feature hits:** the script also requires at least one group's prevalence to be **≥60%**, then computes the larger of positive/negative or negative/positive prevalence ratios, retaining features above `fold_change_threshold`; it saves mean expression among **nonzero** cells as descriptive information.
7. **Search two-feature Boolean gates:** retain features with **>60%** prevalence in at least one group, enumerate feature pairs, and evaluate `AND`, `OR`, `XOR`, `NAND`, `NOR` and `XNOR` using binary expression status (>0). Evaluate each gate separately in receiver-gene-positive and -negative groups, storing gate name, pair, group prevalences and ratio for hits.
8. **Save progressively:** serialize one results dictionary per target receiver gene under `gate_results/` as a pickle; keys are `fold_change` and `gate_analysis`.

## Important code-versus-presentation distinctions

| Topic | Presentation (`logic_search.pdf`) | Current `morphcia.py` |
| --- | --- | --- |
| Spatial representation | Illustrates gridded masks, `grid_size = 5`, 30-unit neighborhoods | Uses a **KD-tree** to select near-donor receivers and direct Euclidean distance around each selected receiver; not the grid-mask helper |
| Organ selector | Slide 24: `check_organ = 'all'` | `check_organ = 'liver'` |
| Parallelism | Slide 24: `use_parallel=True`, `n_jobs=32` | `use_parallel=False`, `n_jobs=2`, `batch_size=100` |
| Neighborhood summary | `max`, `sum` or `mean`; slide defaults to `max` | `cell_environment_representation = 'max'` |
| Minimum prevalence screen | Slide 32: at least 30% in one group | Screen uses **strictly >30%** in one group |
| Boolean search | Discusses one-, two-, and **three-body** interactions | Implements one-feature summaries and **two-feature** gates; no three-feature gate computation |
| Contrast direction | Slide 31 proposes a **symmetric** prevalence-ratio maximum | Two-gate hits currently use **positive/negative** direction (`pos_gate / neg_gate`), not symmetric |
| Statistical tests | Two-proportion z-tests described on slides 9/31 | No statistical p-value calculation is used for the implemented MorphCIA gate hit selection |
| Proposed expression-level rescue | Slide 33 suggests testing mean-expression ratios for prevalent genes | Code reports conditional mean expression but does **not** implement that fallback criterion for single-feature hits |

**Crucial caveat:** In the Boolean pair stage, `(neg_gate == 0)` is marked a hit if prevalence passes the gate minimum, even if positive prevalence is zero. This condition and the one-directional ratio warrant review when interpreting hits; no change to algorithmic behavior has been made in this documentation pass.

## Current configuration and inputs

- `BASE_DIR` points to the original Xenium/mCRC lab filesystem.
- `dist_thresh = 30`, `gene_start_index = 13` and `fold_change_threshold = 1.5`.
- `receiver_cell_type = 'tumor'`; `donor_cell_type = 'fibroblast_WNT5A'`.
- `receiver_gene_names = ['MAP4K4', 'CFB', 'PROX1', 'P4HA1', 'ALDOA', 'TFF3', 'CA12']`.
- Inputs include `sample_organs.txt`, `all_cell_types.pkl`, per-section `morphcia/seurat_to_python/*.csv` discovery paths and corresponding `.feather` tables under a `feathers` directory.
- Uses NumPy, pandas, SciPy, scikit-image, tqdm, joblib and other imports. Several libraries are only imported, not necessarily used.
- The script assumes a fixed **5,083-gene** layout and uses expression columns beginning at position 13, so input schemas must match.

## Interpretation and improvements

Slide 37 shows a proposed AND gate: `fibroblast_WNT5A / C1R` AND `tumor / IQGAP2`, with positive/negative prevalences **0.6885 / 0.4044 (ratio 1.7023)** in colon and **0.6061 / 0.4036 (ratio 1.5016)** in liver. These are presentation examples, **not reproduced results from running the present script**.

Slides 38–39 recommend sweeping gene-expression thresholds (rather than testing only >0), searching distance thresholds, comparing neighborhood aggregation methods and negative-population definitions, constraining rules to include the donor type, and addressing segmentation/spatial-doublet artifacts, temporal discordance, and upregulated versus downregulated targets.

### Running

This is research code with top-level file reads and loops, not a portable CLI. Configure all absolute/relative paths, install dependencies, create the destination output directory and validate cell-type and gene dimensions before running. The PDF and code refer to the same design but **do not represent an identical fixed configuration**.

**Sources:** `logic_search.pdf`, slides 23–39; `morphcia.py`, inspected directly.
