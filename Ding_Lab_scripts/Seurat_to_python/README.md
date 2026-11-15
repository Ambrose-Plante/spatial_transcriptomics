# Seurat to Python data export

## Purpose
Workflow for preparing spatial transcriptomics data from an R/Seurat workflow for downstream Python analysis. The notebook loads a Seurat object saved as a `.qs` file and references a Xenium sample-tracking spreadsheet.

## Inputs and dependencies
- R packages visible in the notebook include `Seurat`, `tidyverse`, `data.table`, `ggplot2`, `gsheet`, `Matrix`, and `qs`.
- A lab-hosted Seurat `.qs` object and an accessible Google Sheet are required.
- Source paths point to the original lab environment and must be adapted.

## Workflow
1. Load Seurat and data-manipulation packages.
2. Read sample metadata from the Xenium tracking sheet.
3. Load the serialized Seurat object.
4. Inspect and prepare spatial cell metadata for downstream processing.
