# Morphological filtering of donor/receiver spatial expression

## Place in the research workflow

`filter_algorithm_morph.ipynb` is an exploratory notebook associated with **Part One** of the accompanying `logic_search.pdf` presentation: *What in a receiver cell is influenced by a donor cell?* (slides 2–22). This README describes the **documented analysis design**, not a verified line-by-line characterization of all code cells in this large notebook.

The goal is to identify receiver-cell genes whose expression prevalence differs between cells **near** a donor-cell population and cells **outside** its neighborhood. The presentation illustrates WNT5A fibroblasts, mCAFs, CD8 tissue-resident-memory T cells, TAMs and tumor cells as donor types, and tumor cells, WNT5A fibroblasts and mCAFs as receiver types.

## Analysis design in the presentation

1. **Rasterize cell locations:** map donor and receiver cell centroids to a shared grid; the illustrated settings are `grid_size = 5` and `distance_threshold = 30`.
2. **Define a near-donor region:** dilate the donor occupancy mask with a disk-shaped neighborhood and intersect with receiver occupancy.
3. **Split the receiver cells/grids:** classify those inside versus outside the dilated donor mask.
4. **Binarize receiver-gene expression:** for every gene in the assay, use 1 if the gene is expressed in a receiver grid and 0 otherwise.
5. **Estimate gene prevalence:** `p_inside = positive_inside / total_inside` and `p_outside = positive_outside / total_outside`.
6. **Compare prevalence:** compute `p_inside / p_outside` and a two-proportion z-test, as outlined on slide 9.
7. **Filter and compare samples:** proposed criteria (slide 10) are `p < 0.05`, ≥10% positivity in either group, and **average fold change ≥1.5** across samples. Slide 11 illustrates per-sample signs (`+`, `-` or `O`), and slides 14–22 show organ-stratified heatmaps.

Example shown on slide 17: proximity to `fibroblast_WNT5A` is associated with higher `MAP4K4` positivity in nearby tumor cells (illustrated ratios ~1.30 colon, 1.25 liver, 1.31 lung, 1.13 breast). These are associations, not evidence of a causal effect.

## Relationship to other files

- [`Seurat_to_python`](../Seurat_to_python/) prepares annotated per-cell expression and coordinates.
- [`fast_morph_scripts`](../fast_morph_scripts/) **does** contain reusable grid mapping, donor-mask dilation and receiver-index extraction functions.
- [`morphcia`](../morphcia/) addresses the **second**, complementary question—candidate nearby cell/gene features associated with a chosen receiver gene's expression.

## Limits and reproducibility

This notebook includes ~6.4 MB of saved notebook material. Its complete cells were **not fully inspected during this README revision**; do **not** assume that each presentation operation above is present or implemented exactly as written in the notebook. In particular, the presentation's slide-specific thresholds are not independently confirmed as notebook defaults.

Review the notebook, adapt original laboratory paths and packages, and verify any produced results against its actual implementation before re-running or citing them. Grid scale, unequal positive/negative group sizes, multiple testing, tissue-specific heterogeneity and segmentation artifacts affect interpretation. Slides 38–39 explicitly flag distance, detection thresholds and spatial doublets as improvement targets.

**Source:** user-provided `logic_search.pdf`, slides 2–22, 38–39; adjacent code in this directory.
