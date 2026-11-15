# Seurat → Python | Spatial-omics data engineering

**Computational bioinformatics · Cross-language interoperability · Xenium / BANKSY / Seurat**

## Why this component exists

Spatial inference is only meaningful if every cell's **identity, tissue section, coordinates, cell-type label, neighborhood annotation, and expression vector** remain aligned. This R-language Jupyter notebook builds that analysis-ready table from a Xenium Seurat object so downstream Python modules can reason over geometry and cell state.

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
