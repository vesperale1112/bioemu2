#!/usr/bin/env python3
"""Evaluate BioEmu2 endpoints on RNAConformer's matched trajectory windows.

Run on an allocated GPU. Output is intentionally kept outside the source tree.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import zipfile
from collections import defaultdict
from io import StringIO
from pathlib import Path

import biotite.structure as bs
import mdtraj as md
import numpy as np
import torch
from biotite.structure.io.pdb import PDBFile
from scipy.spatial import cKDTree

from bioemu2 import EDMSampler, EmuSession, ModifiedResidue, RNA, System, load_model, map_cif_atom_order
from bioemu2.structure.preparation import prepare_sampling_inputs
from bioemu2.trajectory_io import save_mdtraj_cif_frame

RNA_TERMINAL_CCD = {"A": "ADN", "C": "CTN", "G": "GMP", "U": "URI"}


def read_windows(path: Path, limit: int | None) -> list[dict]:
    windows = []
    seen = set()
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            rna = row["rna_id"]
            if rna in seen:
                raise ValueError(f"Only one window per RNA is supported: {rna}")
            seen.add(rna)
            frames = [int(x) for x in row["window_frame_indices"].split(",")]
            if len(frames) < 2 or frames != sorted(set(frames)):
                raise ValueError(f"Invalid window for {rna}: {frames}")
            windows.append({"rna_id": rna, "n_tokens": int(row["n_tokens"]), "frames": frames})
            if limit is not None and len(windows) == limit:
                break
    if not windows:
        raise ValueError("No validation windows selected")
    return windows


def read_manifest(path: Path, rna_ids: set[str]) -> dict[str, dict]:
    result = {}
    with path.open() as stream:
        for line in stream:
            record = json.loads(line)
            rna = record["rna_id"]
            if rna in rna_ids:
                if rna in result:
                    raise ValueError(f"Duplicate manifest entry: {rna}")
                result[rna] = record
    missing = rna_ids - result.keys()
    if missing:
        raise ValueError(f"Missing RNA records: {sorted(missing)}")
    return result


def extract_frames(window: dict, record: dict, archive_root: Path, out: Path) -> list[Path]:
    archive = archive_root / Path(record["archive_path"]).name
    expected_size = record.get("archive_size_bytes")
    if expected_size is not None and archive.stat().st_size != int(expected_size):
        raise ValueError(f"Relocated archive size differs from the original: {archive}")
    references = []
    with zipfile.ZipFile(archive) as source:
        for frame_index in window["frames"]:
            frame = record["frames"][frame_index]
            target = out / "references" / window["rna_id"] / f"frame_{int(frame['time_ps'])}ps.pdb"
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source.read(frame["member_name"]))
            references.append(target)
    return references


def prepare_initial(path: Path, expected_tokens: int):
    atoms = PDBFile.read(StringIO(path.read_text())).get_structure(model=1)
    atoms = atoms[atoms.element != "H"]
    if len(set(atoms.chain_id)) != 1:
        raise ValueError(f"Expected a single RNA chain in {path}")
    atoms.chain_id[:] = "A"
    atoms.atom_name[atoms.atom_name == "O1P"] = "OP1"
    atoms.atom_name[atoms.atom_name == "O2P"] = "OP2"
    starts = bs.get_residue_starts(atoms)
    sequence = "".join(atoms.res_name[starts])
    if len(sequence) != expected_tokens or any(base not in RNA_TERMINAL_CCD for base in sequence):
        raise ValueError(f"Unexpected RNA sequence in {path}: {sequence!r}")
    if len(starts) < 2:
        raise ValueError("At least two residues are required")
    # The original GROMACS RNA has 5'-OH, without P/OP1/OP2. Specify the
    # corresponding unphosphorylated nucleoside so no input atoms are invented.
    first_names = set(atoms.atom_name[:starts[1]])
    if {"P", "OP1", "OP2"} & first_names:
        raise ValueError("Expected an unphosphorylated 5' terminus")
    ccd = RNA_TERMINAL_CCD[sequence[0]]
    atoms.res_name[:starts[1]] = ccd
    system = System([RNA(sequence=sequence, asym_ids=["A"], mods=[ModifiedResidue(res_id=1, ccd_code=ccd)])])
    initial = map_cif_atom_order(system, atoms)
    if initial.n_atoms != len(atoms):
        raise ValueError(f"The mapping dropped atoms: {len(atoms)} input vs {initial.n_atoms} system")
    topology, _ = prepare_sampling_inputs(system, device=torch.device("cpu"))
    if initial.topology != topology:
        raise ValueError("Mapped initial topology differs from the propagation topology")
    has_first_bond = any(
        {a.residue.index, b.residue.index} == {0, 1} and {a.name, b.name} == {"O3'", "P"}
        for a, b in topology.bonds
    )
    if not has_first_bond:
        raise ValueError("First phosphodiester bond is missing")
    if not np.isfinite(initial.xyz).all():
        raise ValueError("Non-finite input coordinates")
    return sequence, ccd, system, initial


def write_scoring_pdb(trajectory: md.Trajectory, sequence: str, output: Path) -> None:
    # The modified terminal CCD is chemically correct; restore the canonical
    # base label only in this scoring copy so USalign recognizes every residue.
    topology = trajectory.topology.copy()
    topology.residue(0).name = sequence[0]
    md.Trajectory(trajectory.xyz.copy(), topology).save_pdb(str(output))


def usalign(binary: Path, prediction: Path, reference: Path) -> dict:
    result = subprocess.run([str(binary), str(prediction), str(reference), "-mol", "RNA", "-outfmt", "2"],
                            capture_output=True, text=True, check=True, timeout=300)
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 5:
            return {"tm_prediction_normalized": float(fields[2]),
                    "tm_reference_normalized": float(fields[3]), "rmsd": float(fields[4]),
                    "aligned_residues": int(fields[10]) if len(fields) > 10 else None}
    raise RuntimeError(f"Cannot parse USalign output: {result.stdout!r} {result.stderr!r}")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def severe_overlap(frame: md.Trajectory) -> dict:
    """Same 1.1 Å non-neighbor-residue diagnostic as the RNAConformer run."""
    xyz = frame.xyz[0] * 10.0
    residues = np.array([atom.residue.index for atom in frame.topology.atoms])
    candidates = cKDTree(xyz).query_pairs(1.1, output_type="ndarray")
    if len(candidates):
        candidates = candidates[np.abs(residues[candidates[:, 0]] - residues[candidates[:, 1]]) > 1]
    clashed = np.unique(candidates) if len(candidates) else []
    return {"severe_nonlocal_overlap_pairs": int(len(candidates)),
            "severe_nonlocal_overlap_atom_fraction": float(len(clashed) / len(xyz))}


def balanced_mean(rows: list[dict], key: str) -> float:
    per_rna = defaultdict(list)
    for row in rows:
        per_rna[row["rna_id"]].append(float(row[key]))
    if not per_rna:
        raise ValueError("Cannot summarize zero rows")
    return float(np.mean([np.mean(values) for values in per_rna.values()]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--usalign", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-windows", type=int)
    parser.add_argument("--max-targets", type=int)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--samples-per-seed", type=int, default=2)
    parser.add_argument("--seeds", type=int, nargs="+", default=[2026, 2027, 2028])
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.max_windows is not None and args.max_windows < 1 or args.max_targets is not None and args.max_targets < 1:
        parser.error("Window/target limits must be positive")
    if min(args.steps, args.samples_per_seed, args.batch_size) < 1 or len(set(args.seeds)) != len(args.seeds):
        parser.error("Steps, sample count and batch size must be positive, and seeds unique")
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"Choose a fresh output directory: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    windows = read_windows(args.windows, args.max_windows)
    manifest = read_manifest(args.manifest, {w["rna_id"] for w in windows})
    prepared = []
    report = {"windows": [], "steps": args.steps, "seeds": args.seeds,
              "samples_per_seed": args.samples_per_seed,
              "note": "Independent endpoints from the same true first frame; short lags are outside the recommended 1–1000 ns.",
              "terminal": "5-prime unphosphorylated nucleoside specified with its CCD; first residue uses modified-residue tokenization."}
    for window in windows:
        record = manifest[window["rna_id"]]
        refs = extract_frames(window, record, args.archive_root, args.output)
        times = [float(record["frames"][i]["time_ps"]) for i in window["frames"]]
        if any(not math.isfinite(t) for t in times) or any(b <= a for a, b in zip(times, times[1:])):
            raise ValueError(f"Non-increasing frame times in {window['rna_id']}: {times}")
        sequence, ccd, system, initial = prepare_initial(refs[0], window["n_tokens"])
        folder = args.output / "inputs" / window["rna_id"]
        folder.mkdir(parents=True)
        save_mdtraj_cif_frame(initial, folder / "initial.cif")
        info = {**window, "sequence": sequence, "times_ps": times, "terminal_ccd": ccd,
                "heavy_atoms": initial.n_atoms}
        report["windows"].append(info)
        prepared.append((window, refs, times, sequence, system, initial))
        print(f"Prepared {window['rna_id']}: {len(sequence)} residues, {initial.n_atoms} atoms, {times} ps", flush=True)
    (args.output / "run_config.json").write_text(json.dumps(report, indent=2) + "\n")
    if args.prepare_only:
        return
    if not torch.cuda.is_available():
        raise RuntimeError("A CUDA GPU is required; submit the supplied PBS job")
    model = load_model(checkpoint_path=args.checkpoint, device="cuda")
    if not model.supports_propagation:
        raise RuntimeError("The selected checkpoint does not support propagation")
    metrics, baselines = [], []
    sampler = EDMSampler(num_steps=args.steps)
    for window, refs, times, sequence, system, initial in prepared:
        rna = window["rna_id"]
        session = EmuSession(system, model=model, mode="equilibrium")
        target_count = min(len(times) - 1, args.max_targets or len(times) - 1)
        for offset in range(1, target_count + 1):
            lag = (times[offset] - times[0]) / 1000.0
            reference = refs[offset]
            baseline = {"rna_id": rna, "frame_index": window["frames"][offset],
                        "time_ps": times[offset], "lag_ns": lag,
                        **usalign(args.usalign, refs[0], reference)}
            baselines.append(baseline)
            write_csv(args.output / "first_frame_baseline.csv", baselines)
            for seed in args.seeds:
                output = args.output / "predictions" / rna / f"seed_{seed}" / f"frame_{window['frames'][offset]}"
                output.mkdir(parents=True)
                generated = session.propagate(initial[[0] * args.samples_per_seed], lag_time_ns=lag,
                                              batch_size=args.batch_size, sampler=sampler,
                                              seed=seed + offset * 100000, progress=False)
                generated.save_xtc(str(output / "samples.xtc"))
                save_mdtraj_cif_frame(generated[0], output / "topology.cif")
                for sample in range(generated.n_frames):
                    prediction = output / f"sample_{sample}.pdb"
                    write_scoring_pdb(generated[sample], sequence, prediction)
                    metrics.append({"rna_id": rna, "frame_index": window["frames"][offset],
                                    "time_ps": times[offset], "lag_ns": lag, "seed": seed,
                                    "sample": sample, "prediction": str(prediction), "reference": str(reference),
                                    **usalign(args.usalign, prediction, reference),
                                    **severe_overlap(generated[sample])})
                write_csv(args.output / "structure_metrics.csv", metrics)
                print(f"Finished {rna} offset={offset} lag_ns={lag:g} seed={seed}", flush=True)
        del session
        torch.cuda.empty_cache()
    summary = {"num_windows": len(prepared), "num_samples": len(metrics),
               "bioemu2": {k: balanced_mean(metrics, k) for k in ("tm_reference_normalized", "rmsd")},
               "first_frame_baseline": {k: balanced_mean(baselines, k) for k in ("tm_reference_normalized", "rmsd")},
               "notes": report["note"]}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
