# RNAConformer matched-window experiment

`examples/evaluate_rna_temporal_bioemu2.py` evaluates BioEmu2 on the windows used by
RNAConformer job `1507016.stdct-mgmt-02`. The 10 IDs and frame indices are read from
RNAConformer's `docs/experiments/flow_structure_20261029/evaluation_windows.csv`, and
the original MD coordinates are read from its trajectory manifest and archives.
Every window starts at 0 ps and has target frames at 100, 200, 300 and 400 ps.

The starting MD frame supplies BioEmu2's full atom coordinates. Each target time is
sampled independently from that same starting frame, not as a time-ordered rollout.
The full run uses 200 EDM steps and three groups of two samples per target; BioEmu2
has a single sampler seed, so these seed groups are not the two independent latent
and decoder noise axes of the RNAConformer experiment. The lags 0.1–0.4 ns lie
below BioEmu2's documented 1–1000 ns recommended range. This is an extrapolation,
not a measurement of BioEmu2 performance in its recommended regime. Conditioning
on all initial-frame coordinates gives BioEmu2 more direct information than the
initial-frame embedding available to RNAConformer.

Original GROMACS RNA uses `O1P/O2P` atom names and an unphosphorylated 5′ end.
The script maps the phosphate names to `OP1/OP2`, defines the first residue as
the matching unphosphorylated nucleoside (CCD `ADN`, `CTN`, `GMP`, or `URI`),
and rejects missing atoms, extra atoms, or a missing first inter-residue bond.
This keeps the complete heavy-atom input without inventing a terminal phosphate.
BioEmu2 tokenizes that first residue as modified; this differs from an ordinary
canonical RNA token and is an additional inference-time representation difference.
The XTC/CIF output retains the terminal CCD. The individual PDB copies used for
USalign relabel only the first residue to the corresponding canonical base, so the
scorer can recognize the full RNA sequence. The coordinates are unchanged.

Per-sample `structure_metrics.csv` uses the same USalign RNA-mode TM/RMSD invocation
and non-neighbor heavy-atom overlap diagnostic as RNAConformer. A separate
`first_frame_baseline.csv` measures the true initial MD structure directly against
each future MD frame; it is not RNAConformer's decoded context baseline. BioEmu2's
built-in physicality check does not cover RNA. The saved coordinates and the
metric tables include every sample; no best-sample selection is applied.

Weights, component files, caches, logs and predictions belong in
`/scratch/csiv212/rnaconformers/bioemu2_rna_20261029`, not the source tree. The
matching `bioemu2-preview` config and checkpoint are expected in `model/`, and
CCD `.cif` files in `ccd/`. After setting these up, submit a short smoke test with
`qsub -v BIOEMU_MODE=smoke scripts/evaluate_rna_temporal_bioemu2.pbs`; the full
job uses `qsub scripts/evaluate_rna_temporal_bioemu2.pbs`. The smoke test uses one
RNA, one target, one sample and 20 diffusion steps, so it is not a benchmark.
Logs are under `runtime/$PBS_JOBID`, and results under `runs/$PBS_JOBID`.
