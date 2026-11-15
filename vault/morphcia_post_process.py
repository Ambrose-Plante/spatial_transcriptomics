import numpy as np
import csv
import time
import os
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import skimage
import glob
import argparse
import gc

# Set up argument parser
parser = argparse.ArgumentParser(description='Morphcia Analysis Script')

# Define command-line arguments
parser.add_argument('--section_id', type=str, required=True, help='Section ID')
parser.add_argument('--gene_grid_filename', type=str, required=True, help='Gene grid filename (relative to base_dir)')
parser.add_argument('--selection_donor', type=str, required=True, help='Donor cell type')
parser.add_argument('--selection_reciever', type=str, required=True, help='Receiver cell type')
parser.add_argument('--base_dir', type=str, required=True, help='Base directory')
parser.add_argument('--grid_size', type=int, default=20, help='Grid size in microns')
parser.add_argument('--max_interaction_dist', type=int, default=100, help='Max interaction distance in microns')
parser.add_argument('--gene_start_index', type=int, default=13, help='Gene start index')
parser.add_argument('--num_cells_expression_filter_cutoff', type=int, default=100, help='Min cells expressing gene')
parser.add_argument('--mean_expression_filter_cutoff', type=float, default=0.001, help='Mean expression filter')
parser.add_argument('--abs_correlation_cutoff', type=float, default=0.2, help='Correlation cutoff to keep pairs')
parser.add_argument('--mode', type=str, default='gene-gene_mask', help='Mode of operation')

# Parse arguments
args = parser.parse_args()

# Ensure required argument section_id is provided
if not args.section_id or not args.selection_donor or not args.selection_reciever or not args.base_dir:
    parser.error('--section_id, --selection_donor, --selection_reciever, and --base_dir are required and must be provided.')

# Assign parsed arguments to variables
section_id = args.section_id
selection_donor = args.selection_donor
selection_reciever = args.selection_reciever
gene_grid_filename = args.gene_grid_filename
BASE_DIR = args.base_dir
grid_size = args.grid_size
max_interaction_dist = args.max_interaction_dist
gene_start_index = args.gene_start_index
num_cells_expression_filter_cutoff = args.num_cells_expression_filter_cutoff
mean_expression_filter_cutoff = args.mean_expression_filter_cutoff
abs_correlation_cutoff = args.abs_correlation_cutoff
mode = args.mode

gene_names_filename = gene_grid_filename.replace("grid", "names")

kernal_size = int(max_interaction_dist/grid_size)

if (kernal_size % 1) != 0:
    raise Exception('max_interaction_dist must be an integer multiple of grid_size')
if (kernal_size % 2) == 0:
    raise Exception('max_interaction_dist must be an odd integer multiple of grid_size')


gene_gridA = np.load(gene_grid_filename)[0]
gene_gridB = np.load(gene_grid_filename)[1]
gene_namesA = np.load(gene_names_filename, allow_pickle=True)
gene_namesB = np.load(gene_names_filename, allow_pickle=True)

gene_inds_keepA = ((gene_gridA.reshape(gene_namesA.shape[0], -1)>0).sum(axis=1))>num_cells_expression_filter_cutoff
gene_gridA = gene_gridA[gene_inds_keepA,:,:]
gene_namesA = gene_namesA[gene_inds_keepA]

gene_inds_keepA = np.mean(gene_gridA.reshape(gene_namesA.shape[0], -1), axis=1)>mean_expression_filter_cutoff
gene_gridA = gene_gridA[gene_inds_keepA,:,:]
gene_namesA = gene_namesA[gene_inds_keepA]

gene_inds_keepB = ((gene_gridB.reshape(gene_namesB.shape[0], -1)>0).sum(axis=1))>num_cells_expression_filter_cutoff
gene_gridB = gene_gridB[gene_inds_keepB,:,:]
gene_namesB = gene_namesB[gene_inds_keepB]

gene_inds_keepB = np.mean(gene_gridB.reshape(gene_namesB.shape[0], -1), axis=1)>mean_expression_filter_cutoff
gene_gridB = gene_gridB[gene_inds_keepB,:,:]
gene_namesB = gene_namesB[gene_inds_keepB]

donor_cell_mask = skimage.morphology.dilation(gene_gridA.sum(axis=0)>0, skimage.morphology.square(kernal_size))
reciever_cell_mask = skimage.morphology.dilation(gene_gridB.sum(axis=0)>0, skimage.morphology.square(kernal_size))

combined_cell_mask = donor_cell_mask * reciever_cell_mask

donor_gene_mask = np.zeros(gene_gridA.shape)
for i in range(donor_gene_mask.shape[0]):
    donor_gene_mask[i,:,:] = skimage.morphology.dilation(gene_gridA[i,:,:], skimage.morphology.square(kernal_size))
donor_gene_mask = donor_gene_mask

reciever_gene_mask = np.zeros(gene_gridB.shape)
for i in range(reciever_gene_mask.shape[0]):
    reciever_gene_mask[i,:,:] = skimage.morphology.dilation(gene_gridB[i,:,:], skimage.morphology.square(kernal_size))
reciever_gene_mask = reciever_gene_mask

Rcell_Dgene_mask = donor_gene_mask * reciever_cell_mask

output_dir = BASE_DIR + "outputs/"+section_id
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

def batch_pearson_correlation(values1, values2):
    n = values1.shape[1]
    mean1 = np.sum(values1, axis=1, keepdims=True) / n
    mean2 = np.sum(values2, axis=1, keepdims=True) / n

    centered1 = values1 - mean1
    centered2 = values2 - mean2

    covariance = np.sum(centered1 * centered2, axis=1)
    variance1 = np.sum(centered1 ** 2, axis=1)
    variance2 = np.sum(centered2 ** 2, axis=1)

    return covariance / (np.sqrt(variance1) * np.sqrt(variance2))

def single_pearson_correlation(values1, values2):
    n = values1.shape[0]
    mean1 = np.sum(values1) / n
    mean2 = np.sum(values2) / n

    centered1 = values1 - mean1
    centered2 = values2 - mean2

    covariance = np.sum(centered1 * centered2)
    variance1 = np.sum(centered1 ** 2)
    variance2 = np.sum(centered2 ** 2)

    return float(covariance / (np.sqrt(variance1) * np.sqrt(variance2)))
    
final_gene_pairs = []

num_donor_genes = gene_gridA.shape[0]
num_reciever_genes = gene_gridB.shape[0]
for k, donor_gene_matrix in enumerate(gene_gridA):
    start_time = time.time()

    if mode=='gene-gene_mask':
        mask = (gene_gridB>0) * donor_gene_mask[k] # grid spaces where there are reciever cells that express geneB, in the proximity (defined by max_interaction_dist) of donor cells that express geneA
        mask = mask>0
        correlations = []
        for rec_gene_ind in range(mask.shape[0]):
            if int((mask[rec_gene_ind]).sum())<6: # less than 6 squares to compare results
                correlations.append(0) # throw out gene pair, not enough information
            else:
                # correlation between the exression of geneB in each receiver cell that expresses it, and the expression of geneA (for geneA>0) in a proximity defined by max_interaction_dist
                correlations.append(single_pearson_correlation(donor_gene_mask[k][mask[rec_gene_ind]], gene_gridB[rec_gene_ind][mask[rec_gene_ind]]))
        correlations = np.asarray(correlations)
        inds_to_keep = np.abs(correlations) > abs_correlation_cutoff
    elif mode=='Dgene-Rcell_mask':
        donor_temp = np.tile(donor_gene_matrix[Rcell_Dgene_mask[k]],(num_reciever_genes,1)) ## repeat donor vector num_reciever_gene times for efficient matrix calculations of correlation
        correlations = batch_pearson_correlation(donor_temp, gene_gridB[:,Rcell_Dgene_mask[k]])
        inds_to_keep = np.abs(correlations) > abs_correlation_cutoff
    else:
        raise Exception('mode must be set')
        
    if inds_to_keep.sum() > 0:
        gene_pair_hits = list(zip(np.repeat(gene_namesA[k], inds_to_keep.sum()).tolist(),
                             gene_namesB[inds_to_keep].tolist(),
                             correlations[inds_to_keep].tolist()))
        #print(gene_pair_hits)
        final_gene_pairs.extend(gene_pair_hits)

    end_time = time.time()
    percent_done = (k + 1) / num_donor_genes * 100
    print(f"compute progress: {percent_done:5.1f} % ({k + 1}/{num_donor_genes}) in {end_time - start_time:5.1f} seconds per iteration")

sorted_gene_pairs = sorted(final_gene_pairs, key=lambda x: abs(x[2]), reverse=True)
#print(sorted_gene_pairs)
df = pd.DataFrame(sorted_gene_pairs, columns=[selection_donor+"Gene", selection_reciever+"Gene", "Correlation"])
df.to_csv(output_dir+"/"+gene_grid_filename.split('/')[-1].split('.npy')[0].split('gene_grid_')[-1] + '.csv', index=False)
