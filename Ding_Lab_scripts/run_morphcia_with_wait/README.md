# Concurrent MorphCIA post-processing launcher

## Purpose in the logic-search research

`run_morphcia_with_wait.sh` is an **execution and scheduling utility**, not the core statistical algorithm. It launches a separate `morphcia_post_process.py` program on numerous **section × donor × receiver** inputs, with limited concurrent subprocesses.

The overarching MorphCIA presentation (`logic_search.pdf`) describes the biological framework: spatial donor/receiver neighborhoods; gene-expression enrichment in close versus distant receiver cells (Part One); and gene/cell-type Boolean-gate features associated with receiver target-gene expression (Part Two). This launcher helps run a *related post-processing workflow* but its named post-processing script is **not included** in this repository, so its detailed scientific output cannot be verified from this Bash file.

## What this script actually does

1. Initialize a Bash shell for Conda and activate the environment named `Morphcia`.
2. Set an explicit list of tissue section IDs (currently only `S16-44514-E7U2`; additional IDs appear in comments).
3. Glob input files at `preprocess/<section_id>/gene_grid*`.
4. Derive donor and receiver labels from the grid filenames and normalize spaces to hyphens.
5. Submit `python morphcia_post_process.py` in the background with identifiers, base directory, gene-grid path, and fixed parameters.
6. Keep at most **10 submitted jobs** in its PID-tracking queue at a time; write each subprocess's standard output/error to a separate `.out` file, then perform a final `wait`.

### Parameters passed to the separate program

| Argument | Value in this script |
| --- | --- |
| `--grid_size` | `20` |
| `--max_interaction_dist` | `100` |
| `--mode` | `gene-gene_mask` |
| `--gene_start_index` | `13` |
| `--num_cells_expression_filter_cutoff` | `100` |
| `--mean_expression_filter_cutoff` | `0.001` |
| `--abs_correlation_cutoff` | `0.2` |

The `20`/`100` geometry settings are **specific to this launcher** and differ from the `5`/`30` settings illustrated in the PDF. They are *not* necessarily equivalent pipelines or parameter scales.

## Requirements and reproducibility

- Bash; a working Conda installation with environment `Morphcia`.
- Prepared `preprocess/<section_id>/gene_grid*` input files.
- The separately referenced `morphcia_post_process.py` and all of its dependencies.
- A valid lab output directory in `--base_dir` (hard-coded to an original personal mount).
- Correct **working directory** for relative paths, especially now that this Bash file resides in its own project folder.

The script limits **concurrent submissions** using PID liveness checks; it does not inspect process exit codes to classify successful versus failed jobs. Verify the `.out` logs independently. Do not assume it executes [`morphcia.py`](../morphcia/morphcia.py): the two scripts are separate.

**Sources:** `run_morphcia_with_wait.sh` and `logic_search.pdf`.
