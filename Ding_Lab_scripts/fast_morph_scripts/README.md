# Fast morphology-based spatial cell utilities

## Purpose
Helper functions for spatial proximity analysis using two-dimensional cell coordinates. The code maps donor and receiver cells onto a shared grid, builds occupancy masks, and identifies receiver cells inside or outside a dilated donor-cell neighborhood.

## Algorithm
- `vectorized_grid_mapping(receiver_cells, donor_cells, grid_size)`: bins `x`/`y` coordinates onto a common grid and returns masks plus mappings back to original cell indices.
- `get_indices_inside_mask(mask, receiver_indices)`: recovers original receiver-cell indices from occupied mask bins.
- `vectorized_split_cells_morph(receiver_cells, donor_cells, grid_size, distance_threshold)`: dilates the donor mask by a disk-shaped structuring element and returns receiver indices classified inside/outside that neighborhood.

## Dependencies and inputs
Requires pandas-like tables with numeric `x` and `y` columns, plus NumPy, pandas, SciPy, scikit-image, matplotlib, statsmodels, and other imported analysis packages. `grid_size` and `distance_threshold` should use consistent coordinate units.

## Usage note
This file defines reusable functions rather than a standalone command-line workflow. Review mask resolution, boundary handling, and the implications of grid-based approximation before interpreting results.
