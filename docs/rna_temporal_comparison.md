# Compare matched RNA validation runs

`examples/compare_rna_temporal_bioemu2.py` compares a completed BioEmu2 run with
RNAConformer job `1507016.stdct-mgmt-02`. It verifies that each seed group has
the same RNA IDs, target frame indices, and samples per target before comparing
RNA-balanced TM-score, USalign RMSD, and severe nonlocal overlap.

Coverage and ensemble precision use the original experiment's 2 Å threshold and
USalign RNA-mode RMSD. Within each RNA and seed group, the 8 generated target
structures are compared with all 4 observed target structures; different seed
groups are never pooled to increase the sample budget. The original RNAConformer
results are read from its existing per-group tables. The BioEmu2 seed groups have
the same sample counts but do not correspond to RNAConformer's two noise axes.
The assessment is descriptive: BioEmu2 is given full initial MD coordinates and
its lag is shorter than its recommended 1–1000 ns range.

Run a CPU comparison after a completed full prediction:

```sh
qsub -v BIOEMU_RUN=/scratch/csiv212/rnaconformers/bioemu2_rna_20261029/runs/JOB_ID \
  scripts/compare_rna_temporal_bioemu2.pbs
```

The output is `comparison/comparison.json` and per-structure coverage/precision
CSVs within the run directory. The PBS job's log is `runtime/JOB_ID/compare.log`.
The comparison exits on missing windows, missing seed groups, or differing sample
counts rather than silently comparing unequal datasets.
