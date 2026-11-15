# MorphCIA spatial neighborhood gene analysis

## Purpose
Analyzes associations between gene expression in receiver cells and gene/cell-type signals in their spatial neighborhood, using Xenium-derived cell tables. The configuration in this script focuses on tumor receiver cells and `fibroblast_WNT5A` donor cells, with organ-specific sample selection.

## Workflow
1. Read a section-to-organ mapping and discover per-section datasets.
2. Select donor and receiver cell populations.
3. Build a `scipy.spatial.cKDTree` to identify receivers near donor cells using a configurable distance threshold.
4. Split nearby receiver cells into positive and negative groups for each receiver gene.
5. Aggregate neighboring cell-type expression profiles (`max`, `sum`, or `mean`).
6. Calculate gene-expression prevalence, fold-change summaries, and candidate gate analyses; save results as serialized data.

## Inputs and dependencies
Uses NumPy, pandas, SciPy, scikit-image, matplotlib, tqdm, joblib, and other Python packages. It expects lab-specific `sample_organs.txt`, `all_cell_types.pkl`, Xenium tables, and hard-coded directory paths.

## Reproducibility
Edit `BASE_DIR`, input paths, cell-type selections, gene list, and thresholds before running. The script executes substantial analysis at top level and assumes the original research data layout; it is not a standalone demo.
