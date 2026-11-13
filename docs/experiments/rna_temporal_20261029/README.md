# BioEmu2 vs RNAConformer: matched RNA windows

This snapshot records BioEmu2 job `1507807.stdct-mgmt-02` and RNAConformer
job `1507016.stdct-mgmt-02`. The 10 RNAs and four target frames per RNA are in
`evaluation_windows.csv` and `bioemu2/run_config.json`. Each method generated
240 structures: three seed groups, two samples per target, with 100, 200, 300
and 400 ps targets. The BioEmu2 samples are independent endpoints from the
true initial structure. RNAConformer's `pure_fm` and `fm_structure` are its
pure flow-matching and structure-fine-tuned models, respectively.

RNA-balanced means from `bioemu2/comparison/comparison.json`:

| Metric | BioEmu2 | RNAConformer pure FM | RNAConformer FM + structure |
| --- | ---: | ---: | ---: |
| Reference-normalized TM-score ↑ | 0.4106 | 0.4428 | 0.4404 |
| USalign RMSD (Å) ↓ | 4.0042 | 3.9028 | 3.9236 |
| Coverage at 2 Å ↑ | 0.0417 | 0.0750 | 0.0750 |
| Ensemble precision at 2 Å ↑ | 0.0625 | 0.1000 | 0.1000 |
| Severe nonlocal overlap atom fraction ↓ | 0.000200 | 0.00000517 | 0.0000103 |

USalign RMSD is the structural-alignment RMSD, not a fixed-residue global
RMSD. Coverage measures the fraction of four reference frames with at least one
of eight predictions within 2 Å; ensemble precision measures the fraction of
predictions within 2 Å of a reference. These are computed per RNA and seed
group, without pooling groups. They cover only the observed window frames.
The overlap diagnostic counts atoms involved in non-neighbor residue heavy-atom
contacts below 1.1 Å; it is not a complete physicality assessment.

The supplementary `bioemu2/comparison/fixed_sequence_c3_rmsd.json` reports
Kabsch RMSD over corresponding C3-prime atoms in sequence order: BioEmu2
6.3101 Å, pure FM 6.1431 Å, FM + structure 6.1915 Å, and the true first-frame
baseline 6.3281 Å. Its source records matching reference-file hashes and a
rigid-transform self-test. Its calculation script is not included in this snapshot.
The USalign first-frame baseline, a separate measure, is in
`bioemu2/first_frame_baseline.csv` and `bioemu2/summary.json`.

These comparisons are descriptive. BioEmu2 receives full initial MD
coordinates while RNAConformer receives initial-frame embeddings. BioEmu2's
0.1–0.4 ns lags are below its recommended 1–1000 ns range, and its single seed
does not correspond to RNAConformer's latent and decoder noise axes. BioEmu2
also treats the unphosphorylated first residue as modified during inference;
see [benchmark details](../../rna_temporal_benchmark.md) and
[comparison method](../../rna_temporal_comparison.md).

The snapshot preserves the original aggregate JSON, per-structure quality,
coverage and precision tables, run configuration, and evaluation windows.
RNAConformer tables include other methods/rows produced by its evaluation;
the comparison selects `temporal_prediction` rows for the matching RNAs.
Absolute scratch paths in the original files record provenance and refer to
structures and other inputs that are not included here. The original structures,
weights and logs remain in scratch; re-running structure alignment requires
those inputs and USalign. Source paths are preserved; CSV line endings were normalized to LF.
