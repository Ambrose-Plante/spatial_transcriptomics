# Fast morphology-based donor/receiver proximity helpers

## Purpose and scientific question

Functions in [`fast_morph_scripts.py`](fast_morph_scripts.py) implement the **spatial masking stage** of *Part One: What in the receiver cell is influenced by the donor cell?* in `logic_search.pdf` (slides 2–10).

The key idea is to compare receiver cells close to a donor-cell type with those outside a predefined donor neighborhood, using gridded occupancy masks and **binary morphological dilation**. This is a *proximity classifier*, not itself an expression-test or causal-inference method.

## How the code works

1. **`vectorized_grid_mapping(receiver_cells, donor_cells, grid_size)`**: builds shared `x_values`/`y_values` grid coordinates from the spatial extent of both groups, bins each cell centroid with `np.searchsorted()`, and creates separate binary donor/receiver occupancy masks. It also retains each receiver's original row index.
2. **`get_indices_inside_mask(mask, receiver_indices)`**: converts occupied mask pixels back into the original receiver table indices.
3. **`vectorized_split_cells_morph(receiver_cells, donor_cells, grid_size, distance_threshold)`**: computes a disk radius of approximately `round(distance_threshold / grid_size)` pixels, dilates the donor mask using `skimage.morphology.binary_dilation`, and intersects it with the receiver mask. Receiver indices falling in the dilated donor mask become `cells_in`; those outside become `cells_out`.

Illustrative presentation settings: `grid_size = 5` and `distance_threshold = 30` (in compatible coordinate units); these are **presentation values**, not hard-coded defaults in these functions.

## Example usage

```python
from fast_morph_scripts import vectorized_split_cells_morph

# donor_cells and receiver_cells are pandas DataFrames with numeric x, y columns
cells_in, cells_out = vectorized_split_cells_morph(
    receiver_cells, donor_cells, grid_size=5, distance_threshold=30
)
near_receivers = receiver_cells.iloc[cells_in]
far_receivers = receiver_cells.iloc[cells_out]
```

The script defines functions only; importing it does not launch a standalone analysis.

## Downstream analysis in the PDF — not implemented in this helper

For each receiver gene, the presentation proposes binary presence/absence scoring among near versus far receiver grids; compute `p_inside` and `p_outside`, fold change `p_inside / p_outside`, and a **two-proportion z-test**; then filter on `p < 0.05`, at least **10%** positivity and average fold change **≥ 1.5** across samples (slides 8–11).

The helper returns cell indexes but **does not compute those expression summaries, tests, filters, or heatmaps**.

## Dependencies and limitations

`numpy`, `pandas`, `scikit-image` and the other libraries imported in the source (including `scipy`, `matplotlib`, `statsmodels`, and `tqdm`) are required for the file as written. Gridding is an approximation: cell collisions in the same bin, grid resolution, edge clipping and disk-radius rounding can affect classifications. Validate coordinate units and results against direct Euclidean-distance calculations if exact-neighbor assignment matters.

**Design reference:** `logic_search.pdf`, slides 3–11. **Implementation reference:** `fast_morph_scripts.py`.
