# Ding Lab spatial transcriptomics: morphological logic search

This directory collects research notebooks and scripts for analyses of spatially resolved gene expression in metastatic colorectal cancer Xenium data. The accompanying design presentation, **`logic_search.pdf` (39 pages; supplied separately, not committed here)**, frames two complementary questions. These are analyses of **spatial association**, not demonstrations of direct causal signaling.

## Analysis 1 — What in a receiver cell is influenced by donor-cell proximity? (slides 2–22)

The presentation compares **receiver** cells located within versus outside a distance threshold of a **donor** cell type.

1. Represent donor and receiver positions on a shared spatial grid. The illustrated settings are `grid_size = 5` and `distance_threshold = 30`.
2. Rasterize donor and receiver occupancy; dilate the donor mask by the distance threshold. Intersect the resulting near/far regions with the receiver mask.
3. For each receiver gene, score expression as present (1) or absent (0) in each selected receiver grid.
4. Compare the fraction of gene-positive receiver grids **inside** and **outside** donor proximity: `p_inside / p_outside`. The presentation uses a two-proportion z-test.
5. Apply the proposed post-processing criteria: `p < 0.05`, at least **10%** gene positivity in either comparison group, and average fold change of at least **1.5** over samples. The figures display results by organ (colon, liver, lung, breast) and paired colon/liver specimens.

Donor-cell types illustrated: `fibroblast_WNT5A`, `mCAF`, `CD8_TRM_T`, `TAM`, and `tumor`. Receiver types illustrated: `tumor`, `fibroblast_WNT5A`, and `mCAF`.

**Related files:** [`fast_morph_scripts`](fast_morph_scripts/) provides the occupancy/morphological-splitting primitives; [`filter_algorithm_morph`](filter_algorithm_morph/) holds an exploratory notebook. The standalone helper script does **not** implement the z-test or sample-level post-processing by itself.

## Analysis 2 — Which nearby gene/cell-type combinations distinguish receiver-gene-positive from negative cells? (slides 23–39)

The presentation reverses the conditioning: start with **receiver** cells positive versus negative for a selected receiver gene (e.g. `PROX1`), and characterize **neighboring cell-type/gene expression**.

1. Select receiver cells near a specified donor population; separate them into **receiver-gene-positive** and **receiver-gene-negative** groups.
2. For each selected receiver cell, summarize expression in surrounding cells across gene × cell-type features (a 3D sample × cell type × gene representation). Example neighborhood aggregation options: `max`, `sum`, `mean`.
3. Compare feature prevalence between positive/negative receiver groups and screen features expressed in at least **30%** of one group.
4. Screen single features and Boolean combinations (e.g. AND, OR, XOR and their complements) for fold-change discrimination. The presentation also *proposes* three-feature Boolean gates.
5. Interpret candidates as **hypotheses** about local microenvironment associations, not established mechanisms.

A displayed example is an **AND** gate combining `fibroblast_WNT5A / C1R` and `tumor / IQGAP2` (ratios about 1.70 in colon and 1.50 in liver; slide 37).

**Related file:** [`morphcia`](morphcia/) implements a specific two-feature Boolean-gate search, with differences from some proposed presentation settings documented in that README.

## Directory guide

| Folder | Role |
| --- | --- |
| [`Seurat_to_python`](Seurat_to_python/) | Export annotated Seurat/Xenium cell data, centroids, and expression into per-section CSV files |
| [`fast_morph_scripts`](fast_morph_scripts/) | Grid occupancy, disk dilation, and near/far cell indexing (Analysis 1 helpers) |
| [`filter_algorithm_morph`](filter_algorithm_morph/) | Historical morphological-filtering notebook (Analysis 1 research context) |
| [`morphcia`](morphcia/) | Receiver-gene-positive/negative neighborhoods and Boolean gate search (Analysis 2) |
| [`run_morphcia_with_wait`](run_morphcia_with_wait/) | Parallel batch launcher for a *separate*, not-included post-processing script |

## Research caveats

- The slides are the **conceptual methods reference**; the repository files are historical implementations, and slide parameters are not automatically the defaults in every script.
- Some code assumes specific lab directory structures, metadata, serialized data, and environments. Inputs and the referenced post-processing script are not all present in this repository.
- Spatial proximity, altered gene prevalence and Boolean gate enrichment do not by themselves establish cell–cell signaling or directionality. Segmentation/spatial doublets, expression detection thresholds, distance choices and temporal mismatch are noted limitations on slides 38–39.
