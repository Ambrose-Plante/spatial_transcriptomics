# MorphCIA | Interpretable Boolean rule mining in spatial tumor microenvironments

**Algorithmic discovery · Structured feature representations · Combinatorial search · Spatial single-cell biology**

## Research question

Can combinations of **neighboring cell types and expressed genes** help explain why nearby receiver cells differ in the expression of a target gene, such as `PROX1`?

[`morphcia.py`](morphcia.py) builds an interpretable spatial feature space and searches for **two-feature Boolean rules** that discriminate **receiver-gene-positive** from **receiver-gene-negative** cells. This is **exploratory rule discovery**, not a trained deep-learning predictor and not evidence of causal intercellular influence.

## Algorithm and representation

**Input:** an annotated spatial cell table per tissue section, with centroids `x,y`, `all_cell_type`, expression columns, and organ/section metadata.

1. **Spatial conditioning with `cKDTree`:** find the nearest donor cell for each receiver; retain receivers within `dist_thresh` of the selected donor type (current donor = `fibroblast_WNT5A`, receiver = tumor classes).
2. **Supervised contrast definition:** for each target receiver gene `t`, label retained receivers `y_i=1` if `expression_i(t)>0`, otherwise `y_i=0`.
3. **Neighborhood representation:** for each selected receiver and each cell-type × gene coordinate, aggregate expression from all cells within the receiver-centered radius using `max`, `sum`, or `mean`. The resulting design array has shape **N receiver cells × C neighborhood cell types × G genes**.
4. **Prevalence-based feature selection:** exclude cell-type/gene combinations absent from both groups at the configured screen (**>30%** positive prevalence in at least one group).
5. **Single-feature association scan:** for features with ≥60% prevalence in at least one group, compute the larger of `P(feature|y=1)/P(feature|y=0)` and its reverse; retain hits above the fold-change threshold. Nonzero conditional mean expression is saved as descriptive metadata.
6. **Pairwise symbolic search:** for feature pairs surviving a stricter **>60%** candidate screen, evaluate **AND, OR, XOR, NAND, NOR, XNOR** on binarized neighborhood expression.
7. **Vectorized scoring and persistence:** evaluate the Boolean gates across positive/negative receiver arrays with NumPy, process pairs in batches (optionally with joblib), and pickle accumulated hits for each target gene.

### Search-space structure

If `K` cell-type/gene features survive screening, the two-feature search covers `K(K−1)/2` unordered pairs × **6** logic gates. Binary gate prevalence is the fraction of receiver neighborhoods satisfying a rule; the script records group prevalences and their **positive/negative ratio** for candidate hits. This explicit, inspectable rule space can be interrogated biologically, unlike opaque latent features.

Current source configuration includes `dist_thresh=30`, `fold_change_threshold=1.5`, `cell_environment_representation='max'`, `check_organ='liver'`, `batch_size=100`, `n_jobs=2`, and `use_parallel=False`. These are **script defaults, not universal analysis settings**.

## Biological example: cross-organ Boolean association

The accompanying methods presentation reports a candidate **AND** rule:

```text
(fibroblast_WNT5A expresses C1R) AND (tumor neighborhood expresses IQGAP2)
```

| Tumor tissue | Target-positive group | Target-negative group | Prevalence ratio |
| --- | ---: | ---: | ---: |
| Colon | 0.6885 | 0.4044 | **1.7023×** |
| Liver | 0.6061 | 0.4036 | **1.5016×** |

These presentation examples illustrate a **potentially conserved, interpretable multicellular neighborhood association**. They are **not fresh results reproduced by executing the present source tree** and should not be described as validated signaling or as predictive-model accuracy.

## Core engineering details

- **Data structures:** `scipy.spatial.cKDTree` for initial donor proximity; NumPy arrays for downstream receiver-centered neighborhood tensors and Boolean evaluation. The later receiver-neighborhood aggregation uses direct distance scans rather than a second KD-tree query.
- **Search optimization:** prevalence gating reduces combinatorial candidates; `itertools.combinations` constructs unordered feature pairs; vectorized Boolean operations process batches of pairs; optional `joblib.Parallel` dispatch is implemented but **disabled by default**.
- **Output contract:** one pickle dictionary under `gate_results/`, with `fold_change` (single-feature descriptors) and `gate_analysis` (gate type, gene/cell-type pair, positive and negative prevalence, ratio) per receiver target gene.
- **Targets in current configuration:** `MAP4K4`, `CFB`, `PROX1`, `P4HA1`, `ALDOA`, `TFF3`, `CA12`. The source assumes **5,083 expression features** and `gene_start_index=13` in its data schema.
