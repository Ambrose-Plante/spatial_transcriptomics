# MorphCIA batch orchestration | Bounded-concurrency research execution

**Scientific workflow automation · Bash · Conda · Batch processing**

## Engineering goal

Spatial transcriptomics studies often require repeated analysis over **tissue sections × donor populations × receiver populations**. [`run_morphcia_with_wait.sh`](run_morphcia_with_wait.sh) is a lightweight **process-scheduling wrapper** that launches a separate MorphCIA post-processing program with bounded local concurrency and per-job logs.

This is an **execution-engineering contribution**, distinct from the spatial-statistics and Boolean-gate algorithms.

## Implemented scheduling pattern

```text
enumerate section IDs
    └─ glob preprocess/<section>/gene_grid*
         └─ parse donor/receiver from filename
              └─ wait for an available PID slot (maximum 10)
                   └─ launch python morphcia_post_process.py ... &
                        └─ redirect stdout/stderr to per-job .out file
final wait for submitted subprocesses
```

The script initializes and activates the `Morphcia` Conda environment; skips nonexistent `gene_grid*` paths; parses and normalizes donor/receiver labels; tracks live PIDs with `kill -0`; and sleeps between queue checks to cap the number of *submitted live processes*.

### Downstream argument contract

| CLI parameter | Script value |
| --- | --- |
| `--grid_size` | `20` |
| `--max_interaction_dist` | `100` |
| `--mode` | `gene-gene_mask` |
| `--gene_start_index` | `13` |
| `--num_cells_expression_filter_cutoff` | `100` |
| `--mean_expression_filter_cutoff` | `0.001` |
| `--abs_correlation_cutoff` | `0.2` |
