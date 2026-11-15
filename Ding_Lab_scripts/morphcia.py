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
from tqdm import tqdm
import pickle
import argparse
from scipy.spatial import cKDTree
from pathlib import Path
from math import comb
from joblib import Parallel, delayed
import multiprocessing

BASE_DIR = '/home/aplante/aplante_diskmnt/mCRC/xenium_5k/Evan/morphcia/'
gene_start_index = 13
dist_thresh = 30 
fold_change_threshold = 1.5
cell_environment_representation = 'max'  # Options: 'max', 'sum', 'mean'
check_organ = 'liver'  # Options: 'colon', 'liver', 'all'
donor_cell_type = 'fibroblast_WNT5A'
receiver_cell_type = 'tumor'
receiver_gene_names = ["MAP4K4", "CFB", "PROX1", "P4HA1", "ALDOA", "TFF3", "CA12"]

use_parallel=False
n_jobs=2
batch_size=100

save_filename = f"{BASE_DIR}gate_results/{donor_cell_type}_{receiver_cell_type}_{check_organ}_{dist_thresh}.pkl"

func_map = {
    'max': np.max,
    'sum': np.sum,
    'mean': np.mean
}
if cell_environment_representation not in func_map:
    raise ValueError(f"Invalid function parameter: {function_name}. Choose from 'max', 'sum', or 'mean'.")
func = func_map[cell_environment_representation]


sample_organs = pd.read_csv('sample_organs.txt', sep='\t')

filenames = glob.glob('morphcia/seurat_to_python/*.csv')
lung_ids = []
colon_ids = []
liver_ids = []
breast_ids = []
section_ids = []
for filename in filenames:
    section_id = filename.split('/')[-1].split('.csv')[0]
    organ = sample_organs[sample_organs['Section ID'] == section_id]['Site']
    if len(organ) != 1:
        raise Exception(f"Unexpected number of matches for {section_id}: {len(organ)}")
    organ = organ.values[0]
    
    if organ == 'Colon':
        colon_ids.append(section_id)
    elif organ == 'Lung':
        lung_ids.append(section_id)
    elif organ == 'Liver':
        liver_ids.append(section_id)
    elif organ == 'Breast':
        breast_ids.append(section_id)
    else:
        print(f"Unknown organ '{organ}' for section_id {section_id}")

    section_ids.append(section_id)
# Optionally, print the collected IDs
print("Colon IDs:", colon_ids)
print("Lung IDs:", lung_ids)
print("Liver IDs:", liver_ids)
print("Breast IDs:", breast_ids)

cross_organ_matches = [
    ["S16-44514-E7U2", "S16-44514-B3U1"],
    ["S15-1909-A5U1", "S15-1909-D2U1"],
    ["S15-32775-D4U2", "S15-32775-A1U2"],
    ["S16-38794-A4U2", "S16-38794-E2U1"],
    ["S13-31378-E7U1", "S13-31378-B1U1"],
    ["HT472C1-S1H3FP1U1", "HT472C1-M1FP1U1"]
]



tumor_cell_types = [
    "APCDD1 tumor",
    "Stem-like tumor",
    "Proliferative stem-like tumor",
    "Proliferative tumor"
]

with open('all_cell_types.pkl', 'rb') as f:
    all_cell_types = pickle.load(f)

# --------------------------------------------------------------------------------------
# This dictionary will store everything:
# {
#   "reciever_gene_1": {
#       "fold_change": [
#           {
#               "gene_name": (cell_type, gene_name),
#               "pos_percent": float,
#               "neg_percent": float,
#               "pos_mean_expr": float or None,
#               "neg_mean_expr": float or None
#           },
#           ...
#       ],
#       "gate_analysis": [
#           {
#               "gate_name": str,
#               "gene1": [cell_type, gene_name],
#               "gene2": [cell_type, gene_name],
#               "pos_val": float,
#               "neg_val": float,
#               "ratio_val": float
#           },
#           ...
#       ]
#   },
#   "reciever_gene_2": {...},
#   ...
# }
# --------------------------------------------------------------------------------------


if check_organ == 'colon':
    ids = colon_ids
elif check_organ == 'liver':
    ids = liver_ids
else:
    ids = section_ids
all_results = {}

for reciever_gene in tqdm(receiver_gene_names):
    print(reciever_gene)
    donor_values_pos = {}
    donor_values_neg = {}

    # ------------------------------------------------------------------------------
    # Process each tissue section
    # ------------------------------------------------------------------------------
    for section_id in tqdm(ids, desc='loading section ids'):
        seurat_table = pd.read_feather('morphcia/seurat_to_python/feathers/' + section_id + '.feather')

        # Partition cells by receiver vs. donor type
        if receiver_cell_type=='tumor':
            receiver_cells = seurat_table[seurat_table["all_cell_type"].isin(tumor_cell_types)]
        else:
            receiver_cells = seurat_table[seurat_table["all_cell_type"] == receiver_cell_type]
        if donor_cell_type=='tumor':
            donor_cells = seurat_table[seurat_table["all_cell_type"].isin(tumor_cell_types)]
        else:
            donor_cells = seurat_table[seurat_table["all_cell_type"] == donor_cell_type]
        
        gene_names = donor_cells.columns[gene_start_index:].values

        if receiver_cells.shape[0] == 0 or donor_cells.shape[0] == 0:
            continue

        # Build KD-Tree for donor coords
        rec_coords = np.stack((receiver_cells["x"].values, receiver_cells["y"].values))
        don_coords = np.stack((donor_cells["x"].values, donor_cells["y"].values))
        don_tree = cKDTree(don_coords.T)

        # Nearest donor neighbor for each receiver cell
        distances, indices = don_tree.query(rec_coords.T, k=1)
        cells_in = np.where(distances < dist_thresh)
        tot_in = cells_in[0].shape[0]
        cells_out = np.where(distances > dist_thresh)
        tot_out = cells_out[0].shape[0]

        if tot_in == 0 or tot_out == 0:
            continue

        rec_near_don = receiver_cells.iloc[cells_in]
        rec_near_don_pos = rec_near_don[rec_near_don[reciever_gene] > 0]
        rec_near_don_neg = rec_near_don[rec_near_don[reciever_gene] == 0]

        rec_coords_pos = np.stack((
            rec_near_don_pos["x"].values,
            rec_near_don_pos["y"].values,
            rec_near_don_pos[reciever_gene].values
        )).T

        rec_coords_neg = np.stack((
            rec_near_don_neg["x"].values,
            rec_near_don_neg["y"].values,
            rec_near_don_neg[reciever_gene].values
        )).T

        # ------------------------------------------------------------------------------
        # POS cells
        # ------------------------------------------------------------------------------
        all_cells_pos = np.zeros((rec_coords_pos.shape[0], all_cell_types.shape[0], gene_names.shape[0]))
        receiver_vals_pos = []
        donor_vals_pos = []

        for i, rec_coord in enumerate(rec_coords_pos):
            x_rec, y_rec = rec_coord[0], rec_coord[1]
            within_thresh = seurat_table[
                np.sqrt((seurat_table['x'] - x_rec) ** 2 + (seurat_table['y'] - y_rec) ** 2) <= dist_thresh
            ]

            for j, cell_type in enumerate(all_cell_types):
                if cell_type=='tumor':
                    temp_cells = within_thresh[within_thresh["all_cell_type"].isin(tumor_cell_types)]
                else:
                    temp_cells = within_thresh[within_thresh["all_cell_type"]==cell_type]
                if len(temp_cells) > 0:
                    all_cells_pos[i, j, :] = func(temp_cells.iloc[:, gene_start_index:].values, axis=0)

            receiver_vals_pos.append(float(rec_coord[2]))

	
        all_cells_pos = np.asarray(all_cells_pos)
        receiver_vals_pos = np.asarray(receiver_vals_pos)

        # ------------------------------------------------------------------------------
        # NEG cells
        # ------------------------------------------------------------------------------
        all_cells_neg = np.zeros((rec_coords_neg.shape[0], all_cell_types.shape[0], gene_names.shape[0]))
        receiver_vals_neg = []
        donor_vals_neg = []

        for i, rec_coord in enumerate(rec_coords_neg):
            x_rec, y_rec = rec_coord[0], rec_coord[1]
            within_thresh = seurat_table[
                np.sqrt((seurat_table['x'] - x_rec) ** 2 + (seurat_table['y'] - y_rec) ** 2) <= dist_thresh
            ]

            for j, cell_type in enumerate(all_cell_types):
                
                if cell_type=='tumor':
                    temp_cells = within_thresh[within_thresh["all_cell_type"].isin(tumor_cell_types)]
                else:
                    temp_cells = within_thresh[within_thresh["all_cell_type"]==cell_type]
                    
                if len(temp_cells) > 0:
                    all_cells_neg[i, j, :] = func(temp_cells.iloc[:, gene_start_index:].values, axis=0)

            receiver_vals_neg.append(float(rec_coord[2]))

        all_cells_neg = np.asarray(all_cells_neg)
        receiver_vals_neg = np.asarray(receiver_vals_neg)

        # ------------------------------------------------------------------------------
        # Check shape and store in dictionaries
        # ------------------------------------------------------------------------------
        if all_cells_pos.shape[2] != 5083 or all_cells_neg.shape[2] != 5083:
            continue

        donor_values_pos[section_id] = all_cells_pos
        donor_values_neg[section_id] = all_cells_neg

    # ------------------------------------------------------------------------------
    print('Concatenate across selected IDs')
    # ------------------------------------------------------------------------------
    if check_organ=='colon':
        array_list = [donor_values_pos[key] for key in donor_values_pos if key in colon_ids]
    elif check_organ=='liver':
        array_list = [donor_values_pos[key] for key in donor_values_pos if key in liver_ids]
    else:
        array_list = [donor_values_pos[key] for key in donor_values_pos]

    positive_temp = np.concatenate(array_list, axis=0) if len(array_list) > 0 else np.array([])

    if check_organ=='colon':
        array_list = [donor_values_neg[key] for key in donor_values_neg if key in colon_ids]
    elif check_organ=='liver':
        array_list = [donor_values_neg[key] for key in donor_values_neg if key in liver_ids]
    else:
        array_list = [donor_values_neg[key] for key in donor_values_neg]

    negative_temp = np.concatenate(array_list, axis=0) if len(array_list) > 0 else np.array([])

    # If either is empty, skip further processing for this gene
    if positive_temp.size == 0 or negative_temp.size == 0:
        all_results[reciever_gene] = {
            "fold_change": [],
            "gate_analysis": []
        }
        continue

    # ------------------------------------------------------------------------------
    print('Filter genes by min expression threshold')
    # ------------------------------------------------------------------------------
    filter_thresh = 30  # minimum percentage expression in either pos/neg
    names = []
    positive_cells = []
    negative_cells = []

    gene_names_df = donor_cells.columns[gene_start_index:].values

    for i, cell_type in enumerate(all_cell_types):
        for j, gene in enumerate(gene_names_df):
            # Skip if it's the receiver gene in receiver cell
            if cell_type == receiver_cell_type and gene == reciever_gene:
                continue

            neg_percent = 100 * (negative_temp[:, i, j] > 0).sum() / negative_temp.shape[0]
            pos_percent = 100 * (positive_temp[:, i, j] > 0).sum() / positive_temp.shape[0]

            if neg_percent > filter_thresh or pos_percent > filter_thresh:
                positive_cells.append(positive_temp[:, i, j])
                negative_cells.append(negative_temp[:, i, j])
                names.append([cell_type, gene])

    positive_cells = np.asarray(positive_cells).T
    negative_cells = np.asarray(negative_cells).T
    names = np.asarray(names)

    if positive_cells.size == 0 or negative_cells.size == 0:
        all_results[reciever_gene] = {
            "fold_change": [],
            "gate_analysis": []
        }
        continue

    neg_percents = 100 * (negative_cells > 0).sum(axis=0) / negative_cells.shape[0]
    pos_percents = 100 * (positive_cells > 0).sum(axis=0) / positive_cells.shape[0]

    # ------------------------------------------------------------------------------
    print('Examine fold change and store results')
    # ------------------------------------------------------------------------------
    fold_change_info = []

    for i in range(len(names)):
        neg_percent = neg_percents[i]
        pos_percent = pos_percents[i]

        # We'll store additional details if there's a fold change
        # or if expression is above a threshold.

        if not (pos_percent < 60 and neg_percent < 60):
            ratio_1 = pos_percent / neg_percent if neg_percent > 0 else np.inf
            ratio_2 = neg_percent / pos_percent if pos_percent > 0 else np.inf
            if np.max((ratio_1, ratio_2)) > fold_change_threshold:
                # Also compute mean expression among those > 0
                pos_mean_expr = 0.0
                neg_mean_expr = 0.0
                pos_nonzero = positive_cells[:, i][positive_cells[:, i] > 0]
                neg_nonzero = negative_cells[:, i][negative_cells[:, i] > 0]

                if pos_nonzero.size > 0:
                    pos_mean_expr = pos_nonzero.mean()
                if neg_nonzero.size > 0:
                    neg_mean_expr = neg_nonzero.mean()

                fold_change_info.append({
                    "gene_name": names[i],
                    "pos_percent": float(round(pos_percent, 2)),
                    "neg_percent": float(round(neg_percent, 2)),
                    "pos_mean_expr": float(round(pos_mean_expr, 2)),
                    "neg_mean_expr": float(round(neg_mean_expr, 2))
                })

    # ------------------------------------------------------------------------------
    print('Identify "consistently expressed" genes')
    # ------------------------------------------------------------------------------
    consistently_expressed_gene_indices = []
    for i in range(len(names)):
        neg_percent = neg_percents[i]
        pos_percent = pos_percents[i]
        if pos_percent > 60 or neg_percent > 60:
            consistently_expressed_gene_indices.append(i)

    consistently_expressed_genes = names[consistently_expressed_gene_indices]

    pos_threshed = positive_cells[:, consistently_expressed_gene_indices] > 0
    neg_threshed = negative_cells[:, consistently_expressed_gene_indices] > 0
    gene_names_subset = np.array(names)[consistently_expressed_gene_indices]

    gate_funcs = {
        "AND":  lambda a, b: a & b,
        "OR":   lambda a, b: a | b,
        "XOR":  lambda a, b: a ^ b,
        "NAND": lambda a, b: ~(a & b),
        "NOR":  lambda a, b: ~(a | b),
        "XNOR": lambda a, b: ~(a ^ b),
    }

    min_fraction = 0.6
    num_genes = pos_threshed.shape[1]
    pairs = np.array(list(itertools.combinations(range(num_genes), 2)))
    total_pairs = pairs.shape[0]

    def process_batch(batch_pairs):
        results = []
        idx1, idx2 = batch_pairs[:, 0], batch_pairs[:, 1]

        pos_1, pos_2 = pos_threshed[:, idx1], pos_threshed[:, idx2]
        neg_1, neg_2 = neg_threshed[:, idx1], neg_threshed[:, idx2]

        for gate_name, gate_func in gate_funcs.items():
            pos_gate = gate_func(pos_1, pos_2).mean(axis=0)
            neg_gate = gate_func(neg_1, neg_2).mean(axis=0)

            valid_mask = (pos_gate >= min_fraction) | (neg_gate >= min_fraction)
            nonzero_mask = neg_gate > 0
            ratio = np.zeros_like(pos_gate)
            ratio[nonzero_mask] = pos_gate[nonzero_mask] / neg_gate[nonzero_mask]

            hits_mask = valid_mask & ((ratio > fold_change_threshold) | (neg_gate == 0))
            hits_idx = np.where(hits_mask)[0]

            for hit in hits_idx:
                i_idx, j_idx = idx1[hit], idx2[hit]
                pos_val = pos_gate[hit]
                neg_val = neg_gate[hit]
                ratio_val = ratio[hit] if neg_gate[hit] > 0 else np.inf

                entry = {
                    "gate_name": gate_name,
                    "gene1": gene_names_subset[i_idx],
                    "gene2": gene_names_subset[j_idx],
                    "pos_val": float(pos_val),
                    "neg_val": float(neg_val),
                    "ratio_val": float(ratio_val),
                }
                results.append(entry)
        return results

    batches = [pairs[start : start + batch_size] for start in range(0, total_pairs, batch_size)]
    if use_parallel:
        parallel_results = Parallel(n_jobs=n_jobs)(
            delayed(process_batch)(batch_pairs) for batch_pairs in tqdm(batches, desc="Processing batches")
        )
    else:
        parallel_results = []
        for batch_pairs in tqdm(batches, desc="Processing batches"):
            parallel_results.append(process_batch(batch_pairs))

    final_results_for_gene = []
    for batch_result in parallel_results:
        final_results_for_gene.extend(batch_result)

    # ------------------------------------------------------------------------------
    print('Store results in our master dictionary under the key for this receiver_gene')
    # ------------------------------------------------------------------------------
    all_results[reciever_gene] = {
        "fold_change": fold_change_info,
        "gate_analysis": final_results_for_gene
    }

    # ------------------------------------------------------------------------------
    print('Save to file (pickle) after each gene iteration')
    # ------------------------------------------------------------------------------
    with open(save_filename, "wb") as f:
        pickle.dump(all_results, f)

