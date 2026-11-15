#!/bin/bash
conda init bash
source ~/.bashrc
conda activate Morphcia

#id_list=( "HT413C1-Th1K2A1U2" "HT472C1-M1FP1U1" "HT472C1-S1H3FP1U1" "S13-18424-B5U1" )

# "S13-31378-B1U1" "S13-31378-E7U1" "S15-1909-A5U1" "S15-1909-D2U1" "S15-21169-A1U1" "S15-23849-A2U1" "S15-32775-A1U2" "S15-32775-D4U2" "S16-38794-A4U2" "S16-38794-E2U1" "S16-44514-B3U1" "S16-44514-E7U2" "S17-16442-AFR1U1" "S20-46186-A16U1" )

id_list=( "S16-44514-E7U2" )

n=10  # Set the max number of parallel jobs
pids=()  # Array to track PIDs of submitted jobs

# Function to wait for at least one job to finish
wait_for_slot() {
    while (( ${#pids[@]} >= n )); do
        for i in "${!pids[@]}"; do
            if ! kill -0 "${pids[i]}" 2>/dev/null; then
                unset 'pids[i]'  # Remove finished job from tracking
            fi
        done
        pids=("${pids[@]}")  # Re-index the array to remove gaps
        sleep 1  # Small delay to avoid excessive CPU usage
    done
}

for section_id in "${id_list[@]}"; do
    for gene_grid_filename in preprocess/"$section_id"/gene_grid*; do
        # Ensure it's a file before proceeding
        [ -f "$gene_grid_filename" ] || continue

        temp1=${gene_grid_filename##*/}
        temp1=${temp1%*.npy}
        reciever=${temp1##*_}
        donor=${temp1%_*}
        donor=${donor##*_}

        # Replace spaces with hyphens
        donor="${donor// /-}"
        reciever="${reciever// /-}"

        echo "submitting: $section_id : $donor : $reciever"

        # Wait if the max job limit is reached
        wait_for_slot

        # Submit the job in the background and capture its PID
        python morphcia_post_process.py \
            --section_id "$section_id" \
            --selection_donor "$donor" \
            --selection_reciever "$reciever" \
            --base_dir "/home/aplante/aplante_diskmnt/mCRC/xenium_5k/Evan/morphcia/" \
            --gene_grid_filename "$gene_grid_filename" \
            --grid_size 20 \
            --max_interaction_dist 100 \
            --mode "gene-gene_mask" \
            --gene_start_index 13 \
            --num_cells_expression_filter_cutoff 100 \
            --mean_expression_filter_cutoff 0.001 \
            --abs_correlation_cutoff 0.2 \
            > "${section_id}_${donor}_${reciever}.out" 2>&1 &

        pid=$!  # Capture PID of the background process
        pids+=("$pid")  # Store the PID in our array
    done
done

# Final wait to ensure all jobs finish before exiting
wait

