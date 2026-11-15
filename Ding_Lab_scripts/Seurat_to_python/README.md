# Seurat-to-Python Xenium cell table export

## Role in the logic-search project

This **R-language Jupyter notebook** prepares annotated single-cell Xenium data for Python-based spatial proximity and Boolean logic analyses described in `logic_search.pdf`. It is upstream **data preparation**, not the statistical logic-search algorithm itself.

The presentation's two analyses require a per-cell spatial location (`x`/`y`), cell type, section identity, and expression values. This notebook assembles those fields from a Seurat object and metadata files.

## Actual notebook workflow

1. Load `Seurat`, `tidyverse`, `data.table`, `ggplot2`, `gsheet`, `Matrix` and `qs`.
2. Read the Xenium tracking Google Sheet and deserialize a lab-hosted `.qs` Seurat object with `qread()`.
3. Read **BANKSY neighborhood/cluster** annotations, align by cell ID, and attach them with `AddMetaData()`.
4. Add **broad cell-type** annotations from a second metadata CSV.
5. Retrieve each section's image-centroid coordinates from `Seurat_object@images[[section_id]]$centroids`.
6. For each section, assemble cell IDs, spatial centroids, metadata, BANKSY cluster and expression values. The expression assay is transposed so **rows are cells and columns are genes**.
7. Write per-section CSV files to `seurat_to_python/<section_id>.csv` using `data.table::fwrite()`; an already existing CSV can be skipped.

The notebook contains another section-export variant later in its cells. Treat the two blocks as historical alternatives, not one validated end-to-end pipeline.

## Outputs and downstream assumptions

Expected CSV fields include cell identifier, spatial centroid coordinates, metadata (including cell-type annotations), and gene-expression columns. Downstream `morphcia.py` assumes per-section tables can be found in its own configured folder, and it uses **Feather** tables in its main analysis; conversion from exported CSVs to that Feather layout is **not documented as implemented in this notebook**.

## Requirements and cautions

- R/IRkernel with the packages above, access to original lab-hosted Xenium data and annotation files, and the tracking sheet.
- Paths in the notebook are lab-specific and need to be replaced before execution.
- The code densifies expression matrices with `as.matrix()`; this may consume substantial memory.
- Check centroid, metadata and expression cell-ID alignment. One current numeric-count check tests nonzero counts rather than enforcing equality, so it should not be relied on as validation.
- Do not share patient-linked metadata or raw clinical sample identifiers without the relevant data governance review.

**Design context:** `logic_search.pdf`, especially slides 3–9 and 24–31; **code inspected:** `Seurat_to_python.ipynb`.
