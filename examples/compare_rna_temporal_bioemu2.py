#!/usr/bin/env python3
"""Compare completed matched-window BioEmu2 and RNAConformer evaluations."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Cannot write an empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def usalign_rmsd(binary: Path, prediction: str, reference: str) -> float:
    result = subprocess.run([str(binary), prediction, reference, "-mol", "RNA", "-outfmt", "2"],
                            capture_output=True, text=True, check=True, timeout=300)
    for line in result.stdout.splitlines():
        if line.strip() and not line.startswith("#"):
            fields = line.split()
            if len(fields) >= 5:
                return float(fields[4])
    raise ValueError(f"Cannot parse USalign output: {result.stdout!r}")


def balanced_summary(rows: list[dict], keys: tuple[str, ...]) -> dict:
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["rna_id"]].append(row)
    if not grouped:
        raise ValueError("Cannot summarize empty results")
    per_rna = {rna: {key: mean(float(row[key]) for row in samples) for key in keys}
               for rna, samples in sorted(grouped.items())}
    return {"rna_count": len(per_rna), "row_count": len(rows),
            "mean": {key: mean(values[key] for values in per_rna.values()) for key in keys},
            "per_rna": per_rna}


def validate_counts(rows: list[dict], expected: set[tuple[str, int]], samples: int, label: str) -> None:
    counts = Counter((row["rna_id"], int(row["frame_index"])) for row in rows)
    if set(counts) != expected or set(counts.values()) != {samples}:
        raise ValueError(f"{label}: target windows or sample counts differ: {counts}")


def original_evaluation(root: Path, method: str, expected: set[tuple[str, int]],
                        samples: int, groups: int, rna_ids: set[str]) -> dict:
    paths = sorted((root / method).glob("latent_*_decoder_*"))
    if len(paths) != groups:
        raise ValueError(f"Expected {groups} seed groups for {method}; found {len(paths)}")
    all_quality, all_coverage, all_precision = [], [], []
    for path in paths:
        quality = [row for row in read_csv(path / "structure_metrics.csv")
                   if row["method"] == "temporal_prediction" and row["rna_id"] in rna_ids]
        validate_counts(quality, expected, samples, f"{method}/{path.name}")
        all_quality.extend(quality)
        coverage = [row for row in read_csv(path / "coverage.csv") if row["rna_id"] in rna_ids]
        precision = [row for row in read_csv(path / "ensemble_precision.csv") if row["rna_id"] in rna_ids]
        # One reference and 'samples' generated frames per target, per seed group.
        expected_refs = Counter(rna for rna, _ in expected)
        if Counter(r["rna_id"] for r in coverage) != expected_refs:
            raise ValueError(f"{method}/{path.name}: coverage references differ")
        expected_preds = Counter({rna: n * samples for rna, n in expected_refs.items()})
        if Counter(r["rna_id"] for r in precision) != expected_preds:
            raise ValueError(f"{method}/{path.name}: precision samples differ")
        all_coverage.extend(coverage)
        all_precision.extend(precision)
    return {
        "quality": balanced_summary(all_quality, ("tm_reference_normalized", "rmsd",
                                                   "severe_nonlocal_overlap_atom_fraction")),
        "coverage": balanced_summary(all_coverage, ("covered", "nearest_rmsd")),
        "precision": balanced_summary(all_precision, ("within_reference_set", "nearest_rmsd")),
    }


def calculate_ensembles(rows: list[dict], binary: Path, threshold: float) -> tuple[list[dict], list[dict]]:
    windows = defaultdict(list)
    for row in rows:
        windows[row["rna_id"], row["seed"]].append(row)
    coverage, precision = [], []
    for (rna, seed), group in sorted(windows.items()):
        predictions = sorted({row["prediction"] for row in group})
        references = sorted({row["reference"] for row in group})
        if len(predictions) != len(group):
            raise ValueError(f"Repeated prediction paths for {rna}, seed {seed}")
        distances = {(p, r): usalign_rmsd(binary, p, r)
                     for p, r in itertools.product(predictions, references)}
        for reference in references:
            distance = min(distances[p, reference] for p in predictions)
            coverage.append({"rna_id": rna, "seed": seed, "reference": reference,
                             "nearest_rmsd": distance, "covered": float(distance <= threshold)})
        for prediction in predictions:
            distance = min(distances[prediction, r] for r in references)
            precision.append({"rna_id": rna, "seed": seed, "prediction": prediction,
                              "nearest_rmsd": distance,
                              "within_reference_set": float(distance <= threshold)})
        print(f"Compared ensemble {rna}, seed {seed}: {len(predictions)} predictions, {len(references)} references", flush=True)
    return coverage, precision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bioemu-results", type=Path, required=True)
    parser.add_argument("--rna-results", type=Path, required=True)
    parser.add_argument("--usalign", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads((args.bioemu_results / "run_config.json").read_text())
    original = json.loads((args.rna_results / "comparison.json").read_text())
    samples = int(config["samples_per_seed"])
    seeds = config["seeds"]
    if samples != original["config"]["latent_samples"] * original["config"]["decoder_samples"]:
        raise ValueError("Samples per seed differ from the RNAConformer evaluation")
    expected = {(w["rna_id"], int(frame)) for w in config["windows"] for frame in w["frames"][1:]}
    rna_ids = {rna for rna, _ in expected}
    rows = read_csv(args.bioemu_results / "structure_metrics.csv")
    for seed in seeds:
        validate_counts([r for r in rows if int(r["seed"]) == seed], expected, samples, f"BioEmu2 seed {seed}")
    if len(rows) != len(expected) * len(seeds) * samples:
        raise ValueError("BioEmu2 produced unexpected extra rows")
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"Choose a fresh output directory: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    threshold = float(original["config"]["coverage_rmsd_angstrom"])
    coverage, precision = calculate_ensembles(rows, args.usalign, threshold)
    write_csv(args.output / "bioemu2_coverage.csv", coverage)
    write_csv(args.output / "bioemu2_ensemble_precision.csv", precision)
    results = {
        "bioemu2": {
            "quality": balanced_summary(rows, ("tm_reference_normalized", "rmsd",
                                                "severe_nonlocal_overlap_atom_fraction")),
            "coverage": balanced_summary(coverage, ("covered", "nearest_rmsd")),
            "precision": balanced_summary(precision, ("within_reference_set", "nearest_rmsd")),
        },
        "pure_fm": original_evaluation(args.rna_results, "pure_fm", expected, samples, len(seeds), rna_ids),
        "fm_structure": original_evaluation(args.rna_results, "fm_structure", expected, samples, len(seeds), rna_ids),
    }
    comparison = {
        "rmsd_threshold_angstrom": threshold,
        "rna_results": str(args.rna_results),
        "bioemu_results": str(args.bioemu_results),
        "methods": results,
        "notes": ["Each seed group has the same number of samples per window; groups are not pooled for coverage.",
                  "BioEmu2 uses actual first-frame coordinates and samples independent endpoints; RNAConformer uses first-frame embeddings.",
                  "The 0.1–0.4 ns BioEmu2 lags are below its recommended training range.",
                  "USalign RMSD is the alignment RMSD, and coverage is limited to these observed window frames."],
    }
    (args.output / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    print(json.dumps({method: {name: metric["mean"] for name, metric in values.items()}
                      for method, values in results.items()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
