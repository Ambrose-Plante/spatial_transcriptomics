# Seurat → Python | Spatial-omics data engineering

**Computational bioinformatics · Cross-language interoperability · Xenium / BANKSY / Seurat**

## Why this component exists

Spatial inference is only meaningful if every cell's **identity, tissue section, coordinates, cell-type label, neighborhood annotation, and expression vector** remain aligned. This R-language Jupyter notebook builds that analysis-ready table from a Xenium Seurat object so downstream Python modules can reason over geometry and cell state.

This is the **ETL and data-contract layer** of the [Ding Lab spatial logic-search project](../README.md), not the statistical test or Boolean classifier itself.

## Computational workflow

1. **Load annotated single-cell data:** deserialize a `.qs` Seurat object and read metadata/section-tracking inputs.
2. **Join neighborhood information:** align external **BANKSY** neighborhood/cluster annotations to Seurat cell identifiers and attach with `AddMetaData()`.
3. **Resolve biological classes:** join broad cell-type labels for donor/receiver stratification.
4. **Recover geometry:** extract per-section spatial centroids from `Seurat_object@images[[section_id]]$centroids`.
5. **Materialize expression features:** extract the selected assay's expression matrix; transpose from **genes × cells** to **cells × genes**.
6. **Export section-level tables:** combine `cell_id`, `x`/`y` centroids, section metadata, cell-type/BANKSY annotations, and gene columns; write `seurat_to_python/<section>.csv` via `data.table::fwrite()`.

### Downstream data contract

```text
one row = one cell
[cell_id | x | y | section/annotation columns | gene_1 | ... | gene_G]
```

The MorphCIA scripts consume annotated spatial cell tables with specific column offsets and category names. **This notebook exports CSVs; it does not implement the later CSV→Feather conversion** assumed by `morphcia.py`.

## Engineering decisions & trade-offs

- **Cell-ID joins are the invariant:** coordinate rows, metadata rows and assay columns must refer to the same cells in the same order. Rownames, section IDs and `orig.ident` need explicit validation on each dataset.
- **Per-section partitioning** reduces bookkeeping complexity in multi-tissue studies and permits independent organ-level analyses.
- **Dense conversion is expensive:** `as.matrix()` can expand sparse transcript counts into a large in-memory matrix. For larger cohorts, chunked or sparse-friendly Arrow/Parquet exports would be advisable; **these optimizations are recommendations, not features claimed here**.
- The notebook contains **two historical export variants**; they are not a versioned, productionized pipeline. Its current numeric-count check does not reliably enforce equality and should be replaced with an explicit assertion before reuse.

## What this demonstrates

Reproducible handling of single-cell identifiers and annotations; R/Seurat-to-Python interoperation; matrix reorientation; spatial coordinate extraction; and research data engineering for downstream high-dimensional algorithms.

**Source:** [`Seurat_to_python.ipynb`](Seurat_to_python.ipynb). **Dependencies:** R/IRkernel, Seurat, tidyverse, data.table, Matrix, qs, gsheet and lab-hosted source data. Adapt hard-coded paths and respect clinical-data access restrictions.
