#!/bin/bash
conda init bash
source ~/.bashrc
conda activate Morphcia


#id_list=( "HT413C1-Th1K2A1U2" "HT472C1-M1FP1U1" "HT472C1-S1H3FP1U1" "S13-18424-B5U1" )

# "S13-31378-B1U1" "S13-31378-E7U1" "S15-1909-A5U1" "S15-1909-D2U1" "S15-21169-A1U1" "S15-23849-A2U1" "S15-32775-A1U2" "S15-32775-D4U2" "S16-38794-A4U2" "S16-38794-E2U1" "S16-44514-B3U1" "S16-44514-E7U2" "S17-16442-AFR1U1" "S20-46186-A16U1" )

id_list=( "S16-44514-E7U2" )

for section_id in "${id_list[@]}"; do
	for gene_grid_filename in `ls preprocess/${section_id}/gene_grid*`; do
		temp1=${gene_grid_filename##*/}
		temp1=${temp1%*.npy}
		reciever=${temp1##*_}
		donor=${temp1%_*}
		donor=${donor##*_}
		echo "submitting: $section_id : $donor : $reciever"
		python morphcia_post_process.py --section_id $section_id --selection_donor $donor --selection_reciever $reciever --base_dir "/home/aplante/aplante_diskmnt/mCRC/xenium_5k/Evan/morphcia/" --gene_grid_filename $gene_grid_filename --grid_size 20 --max_interaction_dist 100 --mode "gene-gene_mask" --gene_start_index 13 --num_cells_expression_filter_cutoff 5 --mean_expression_filter_cutoff 0.001 --abs_correlation_cutoff 0.2 > ${section_id}_${donor}_${reciever}.out 2>&1 &

	done

done


