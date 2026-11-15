# Parallel MorphCIA post-processing launcher

## Purpose
Bash launcher for running multiple spatial gene-interaction post-processing jobs while limiting the number of simultaneous processes.

## Workflow
1. Initialize Conda and activate the `Morphcia` environment.
2. Iterate selected tissue section IDs and `preprocess/<section_id>/gene_grid*` files.
3. Parse donor and receiver identifiers from each filename.
4. Launch `morphcia_post_process.py` in the background with section, grid, distance, and expression-filter parameters.
5. Track process IDs, limit concurrency to `n=10`, and wait for all jobs to finish. Each job writes a separate `.out` log.

## Requirements
Bash, Conda environment `Morphcia`, prepared gene-grid files, and the separate `morphcia_post_process.py` script (not included in this folder). The launcher includes a lab-specific absolute base path.

## Reproducibility
Update the section list, base directory, environment name, and referenced post-processing script before use. Because this file is now inside its own subfolder, run it from the intended project working directory or adjust relative paths.
