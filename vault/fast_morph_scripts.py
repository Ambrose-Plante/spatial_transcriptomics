import numpy as np
import csv
import time
import os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import skimage
import glob
from matplotlib.lines import Line2D
import tqdm
import pickle
from scipy.spatial import cKDTree
from statsmodels.stats.proportion import proportions_ztest
from skimage.morphology import binary_dilation, disk

def vectorized_grid_mapping(
    receiver_cells, 
    donor_cells, 
    grid_size
):
    # Compute grid extents
    x_min = np.floor(min(donor_cells["x"].min(), receiver_cells["x"].min()))
    x_max = np.ceil(max(donor_cells["x"].max(), receiver_cells["x"].max()))
    y_min = np.floor(min(donor_cells["y"].min(), receiver_cells["y"].min()))
    y_max = np.ceil(max(donor_cells["y"].max(), receiver_cells["y"].max()))
    
    # Build x and y "edges" for the bins
    x_values = np.arange(x_min - grid_size, x_max + grid_size, grid_size)
    y_values = np.arange(y_min - grid_size, y_max + grid_size, grid_size)
    
    # We assume x_values is sorted ascending, likewise y_values
    nx = len(x_values)
    ny = len(y_values)
    
    # Extract coordinates
    rec_coords = np.column_stack([receiver_cells["x"].values, receiver_cells["y"].values])
    don_coords = np.column_stack([donor_cells["x"].values, donor_cells["y"].values])
    
    # Initialize 2D masks
    donor_cell_mask = np.zeros((nx, ny), dtype=int)
    receiver_cell_mask = np.zeros((nx, ny), dtype=int)

    # -----------------------------------------------------------------------
    # 1) Map each receiver coordinate to the proper bin (x_idx, y_idx)
    #    using np.searchsorted, side='right' - 1 => bin index
    # -----------------------------------------------------------------------
    rx = np.searchsorted(x_values, rec_coords[:, 0], side='right') - 1
    ry = np.searchsorted(y_values, rec_coords[:, 1], side='right') - 1

    # Fix potential out-of-bounds (e.g., coordinate < x_values[0])
    rx = np.clip(rx, 0, nx - 1)
    ry = np.clip(ry, 0, ny - 1)

    # For the mask: set each corresponding bin to 1
    #   (If multiple points map to the same bin, it remains 1.)
    receiver_cell_mask[rx, ry] = 1

    # Build receiver_indices as an array of [x_idx, y_idx, original_index]
    # np.arange(len(rx)) gives the original cell's index
    receiver_indices = np.column_stack([rx, ry, np.arange(len(rx))])

    # -----------------------------------------------------------------------
    # 2) Repeat for donor coordinates
    # -----------------------------------------------------------------------
    dx = np.searchsorted(x_values, don_coords[:, 0], side='right') - 1
    dy = np.searchsorted(y_values, don_coords[:, 1], side='right') - 1

    dx = np.clip(dx, 0, nx - 1)
    dy = np.clip(dy, 0, ny - 1)

    donor_cell_mask[dx, dy] = 1
    donor_indices = np.column_stack([dx, dy, np.arange(len(dx))])

    return (
        x_values, y_values,
        receiver_cell_mask,
        donor_cell_mask,
        receiver_indices,
        donor_indices
    )


def get_indices_inside_mask(mask, receiver_indices):
    # 1) Identify (i, j) where both masks are nonzero
    row_idxs, col_idxs = np.where(mask > 0)

    # 2) Build a dictionary: (x_idx, y_idx) -> list of receiver cell indices
    bin_map = {}
    for x_idx, y_idx, rec_idx in receiver_indices:
        bin_map.setdefault((x_idx, y_idx), []).append(rec_idx)

    # 3) Collect all matching receiver indices
    cells = []
    for i, j in zip(row_idxs, col_idxs):
        cells.extend(bin_map.get((i, j), []))

    return cells


def vectorized_split_cells_morph(receiver_cells, donor_cells, grid_size, distance_threshold):
	x_values,y_values,receiver_cell_mask,donor_cell_mask,receiver_indices,donor_indices=vectorized_grid_mapping(receiver_cells,donor_cells,grid_size)
	radius_pixels = int(round(distance_threshold / grid_size))
	selem = disk(radius_pixels)
	donor_cell_dilated_mask = binary_dilation(donor_cell_mask, selem)
	inverted_donor_cell_dilated_mask = 1-donor_cell_dilated_mask
	
	cells_in = get_indices_inside_mask(donor_cell_dilated_mask * receiver_cell_mask, receiver_indices)
	cells_out = get_indices_inside_mask(inverted_donor_cell_dilated_mask * receiver_cell_mask, receiver_indices)

	return cells_in, cells_out	
