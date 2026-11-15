# Fast spatial morphology | Donor–receiver proximity indexing

**Algorithm design · Geometric feature engineering · Vectorized scientific Python**

## Research problem

How can millions of spatial **cell–gene observations** be stratified efficiently by their location relative to a biologically defined donor population?

[`fast_morph_scripts.py`](fast_morph_scripts.py) supplies the geometric building blocks for **near-versus-far receiver-cell comparisons**. Rather than computing a full donor-by-receiver pairwise-distance matrix, it rasterizes cell centroids and uses **binary morphological dilation** to build a discretized donor-neighborhood mask.

## Implemented functions

| Function | Algorithm / output |
| --- | --- |
| `vectorized_grid_mapping(receiver_cells, donor_cells, grid_size)` | Determines shared spatial extents; bins `x, y` using `np.searchsorted`; returns donor/receiver occupancy matrices and per-cell bin/row-index mappings |
| `get_indices_inside_mask(mask, receiver_indices)` | Maps occupied receiver bins back to **all original receiver rows**, preserving multiple cells that fall in the same grid bin |
| `vectorized_split_cells_morph(receiver_cells, donor_cells, grid_size, distance_threshold)` | Dilates donor occupancy by a disk with pixel radius `round(distance_threshold / grid_size)`, intersects with receiver occupancy, and returns **near/far row-index lists** |

### Algorithmic idea

```text
donor centroids  ──► 2-D occupancy grid ──► disk dilation ──┐
                                                             ├─► receiver indices: near / far
receiver centroids ─► 2-D occupancy grid ──► mask intersect ──┘
```

This is a **grid-based approximation to a distance query**. By representing a spatial neighborhood as a binary mask, selection becomes image morphology plus index lookup. Memory depends on grid extent/resolution; fine grids, crowded bins and radius rounding affect geometry and accuracy.

## Python example

```python
from fast_morph_scripts import vectorized_split_cells_morph

# pandas DataFrames with numeric x and y columns in the same coordinate frame
near_idx, far_idx = vectorized_split_cells_morph(
    receiver_cells, donor_cells, grid_size=5, distance_threshold=30
)
near_receivers = receiver_cells.iloc[near_idx]
far_receivers  = receiver_cells.iloc[far_idx]
```

The `5`/`30` values are **illustrative research-presentation settings**, not hard-coded library defaults. Validate coordinate units and occupancy behavior against exact Euclidean or KD-tree queries when precision is important.

## Why this matters biologically

The near/far partition defines a spatially meaningful contrast for testing whether a **receiver gene's detection prevalence** differs in a donor cell type's local microenvironment. The accompanying [filtering notebook](../filter_algorithm_morph/) explores organ-stratified differences in metastatic colorectal cancer Xenium sections. See the [WNT5A-fibroblast neighborhood result](../README.md#biological-results-organ-stratified-spatial-associations).

**Boundary of responsibility:** this file implements geometry and cell-index recovery. It imports `proportions_ztest` but **does not run a gene-level z-test, calculate near/far fold changes, apply significance filters, or generate biological heatmaps**. Those belong to the broader exploratory design, not to the three helper functions.
