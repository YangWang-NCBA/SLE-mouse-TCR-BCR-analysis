#!/usr/bin/env python3
"""Process paired-end TCR/BCR libraries with Trim Galore and MiXCR 4.7.0."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from sequence_metrics import assemble_chain_counts

CODE = Path(__file__).resolve().parents[1]
TAG = (r"^N{0:2}aagcagtggtatcaacgcagagt(UMI:N{4}TN{4}TN{4})tcttgggg(R1:*) "
       r"\ ^(R2:*) || ^(R1:*) \ "
       r"^N{0:2}aagcagtggtatcaacgcagagt(UMI:N{4}TN{4}TN{4})tcttgggg(R2:*)")
CHAINS = {"BCR": ("IGH", "IGK", "IGL"), "TCR": ("TRA", "TRB", "TRG", "TRD")}
FEATURES = ("VRegionTrimmed", "FR1", "CDR1", "FR2", "CDR2", "FR3")
def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_inputs(manifest, raw_root, selected=None):
    """Locate the R1 and R2 files for each selected library."""
    rows = read_tsv(manifest)
    if not rows:
        raise ValueError("The library manifest has no libraries.")
    needed = {"analysis_library_id", "receptor", "R1_filename", "R2_filename"}
    if not needed <= rows[0].keys():
        raise ValueError("The manifest requires: " + ", ".join(sorted(needed)))
    identities, names = set(), set()
    for row in rows:
        ident = row["analysis_library_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", ident) or ident in identities:
            raise ValueError("Unsafe or duplicate analysis_library_id.")
        identities.add(ident)
        if row["receptor"] not in CHAINS:
            raise ValueError("Unknown receptor category for " + ident)
        for direction in ("R1", "R2"):
            name = row[direction + "_filename"]
            if not name or "/" in name or "\\" in name or name in names:
                raise ValueError("Unsafe or repeated FASTQ basename: " + name)
            names.add(name)
    if selected:
        if selected not in identities:
            raise ValueError("Library ID is not in the manifest: " + selected)
        rows = [r for r in rows if r["analysis_library_id"] == selected]
    chosen = {row[d + "_filename"] for row in rows for d in ("R1", "R2")}
    found = {name: [] for name in chosen}
    for path in Path(raw_root).rglob("*"):
        if path.is_file() and path.name in found:
            found[path.name].append(path.resolve())
    inputs = []
    for row in rows:
        mates = []
        for direction in ("R1", "R2"):
            name = row[direction + "_filename"]
            if len(found[name]) != 1 or found[name][0].stat().st_size == 0:
                raise ValueError("Expected one nonempty FASTQ: " + name)
            mates.append(found[name][0])
        inputs.append((row, mates))
    return inputs


def commands(row, mates, root, mixcr, trim, threads):
    ident = row["analysis_library_id"]
    clean = root / "clean" / ident
    folder = root / "results" / ident
    prefix = folder / ident
    aligned, refined, clns = [Path(str(prefix) + ext) for ext in (".vdjca", ".refined.vdjca", ".clns")]
    prepared = [clean / (ident + "_val_" + str(n) + ".fq.gz") for n in (1, 2)]
    stage = [
        ("trim", [trim, "--paired", "--cores", "1", "--quality", "20", "--length", "35",
                  "--fastqc", "--gzip", "--basename", ident, "--output_dir", str(clean),
                  str(mates[0]), str(mates[1])]),
        ("align", [mixcr, "align", "-t", str(threads), "--species", "mmu", "--library", "imgt.202312-3.sv8",
                   "--rna", "-p", "generic-amplicon-with-umi", "--assemble-clonotypes-by", "VDJRegion",
                   "-OvParameters.geneFeatureToAlign=VRegionWithP", "--tag-pattern", TAG,
                   "--not-aligned-R1", str(prefix) + "_1_val_1un.fq.gz",
                   "--not-aligned-R2", str(prefix) + "_2_val_2un.fq.gz",
                   "--floating-left-alignment-boundary", "--floating-right-alignment-boundary", "C",
                   "--json-report", str(prefix) + "_align.report.json", "--report", str(prefix) + "_align.report.txt",
                   str(prepared[0]), str(prepared[1]), str(aligned)]),
        ("refine", [mixcr, "refineTagsAndSort", "--report", str(prefix) + ".refine.report.txt", str(aligned), str(refined)]),
        ("assemble", [mixcr, "assemble", "--assemble-clonotypes-by", "VDJRegion",
                      "-OseparateByJ=true", "-OseparateByV=true", "-OseparateByC=true",
                      "--report", str(prefix) + ".assemble.report.txt", str(refined), str(clns)]),
    ]
    if row["receptor"] == "BCR":
        stage.append(("findAlleles", [mixcr, "findAlleles", "-t", str(threads), "--verbose",
                      "--report", str(prefix) + ".findAlleles.report.txt",
                      "--export-alleles-mutations", str(prefix) + ".alleles.tsv",
                      "--export-library", str(prefix) + ".alleles.json",
                      "--output-template", str(folder / "{file_name}.allelic.clns"), str(clns)]))
        clns = Path(str(prefix) + ".allelic.clns")
    for chain in CHAINS[row["receptor"]]:
        cmd = [mixcr, "exportClones", "--chains", chain, "--not-covered-as-empty", "-chains",
               "-allNFeatures", "FR1Begin", "FR4End", "-allAAFeatures", "FR1Begin", "FR4End"]
        if row["receptor"] == "BCR":
            if chain == "IGH":
                cmd += ["-isotype", "auto"]
            for feature in FEATURES:
                cmd += ["-nMutationsCount", feature, "substitutions", "-nMutationsRate", feature, "substitutions"]
        cmd += [str(clns), str(prefix) + "_" + chain + ".tsv"]
        stage.append(("export_" + chain, cmd))
    return stage


def execute_library(item, args, root, mixcr, trim):
    row, mates = item
    ident = row["analysis_library_id"]
    log = root / "work" / ident
    log.mkdir(parents=True)
    (root / "clean" / ident).mkdir(parents=True)
    (root / "results" / ident).mkdir(parents=True)
    plan = commands(row, mates, root, mixcr, trim, args.threads)
    environment = os.environ.copy()
    environment["JAVA_TOOL_OPTIONS"] = (environment.get("JAVA_TOOL_OPTIONS", "") +
                                        " -XX:ActiveProcessorCount=" + str(args.threads)).strip()
    for stage, argv in plan:
        if stage.startswith("export_"):
            report_path = root / "results" / ident / (ident + ".assemble.report.txt")
            report_text = report_path.read_text(encoding="utf-8")
            counts, closed = assemble_chain_counts(report_text)
            if closed and counts[stage[7:]][0] == 0:
                continue
        with (log / (stage + ".log")).open("w", encoding="utf-8") as stream:
            completed = subprocess.run(argv, stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
                                       env=environment, check=False)
        if completed.returncode:
            raise RuntimeError(f"{ident}: {stage} failed. See work/{ident}/{stage}.log")
        if stage == "trim":
            for n in (1, 2):
                prepared = root / "clean" / ident / (ident + "_val_" + str(n) + ".fq.gz")
                if not prepared.is_file() or not prepared.stat().st_size:
                    raise RuntimeError(ident + ": expected trimmed mate was not produced.")
        if stage.startswith("export_"):
            exported = Path(argv[-1])
            if not exported.is_file() or not exported.stat().st_size:
                raise RuntimeError(f"MiXCR did not produce {exported.name}")


def find_tool(value):
    resolved = shutil.which(value)
    if not resolved:
        raise ValueError("Required executable not found: " + value)
    return resolved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path, help="New processing directory; must not exist")
    parser.add_argument("--manifest", type=Path, required=True, help="Prepared data/metadata/libraries.tsv")
    parser.add_argument("--mixcr", default="mixcr", help="MiXCR 4.7.0 executable")
    parser.add_argument("--trim-galore", default="trim_galore")
    parser.add_argument("--threads", type=int, default=16, help="MiXCR threads per library")
    parser.add_argument("--jobs", type=int, default=1, help="Concurrent libraries")
    parser.add_argument("--library-id", help="Process one manifest-selected library")
    args = parser.parse_args()
    try:
        root, raw = args.out_root.resolve(), args.raw_root.resolve()
        if root.exists() or root == raw or raw in root.parents or root in raw.parents or root == CODE or CODE in root.parents or root in CODE.parents:
            raise ValueError("Use a new output directory outside the raw-data and code directories.")
        if not raw.is_dir():
            raise ValueError("The raw FASTQ directory does not exist.")
        if args.threads < 8 or args.jobs < 1:
            raise ValueError("Use at least 8 threads per library and at least one job (paired trimming overhead).")
        allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
        if allocated and args.jobs * args.threads > allocated:
            raise ValueError("jobs times threads exceeds the Slurm CPU allocation.")
        items = load_inputs(args.manifest, raw, args.library_id)
        mixcr, trim = find_tool(args.mixcr), find_tool(args.trim_galore)
        version = subprocess.run([mixcr, "--version"], capture_output=True, text=True, check=True)
        if not re.search(r"\b4\.7\.0\b", version.stdout + version.stderr):
            raise ValueError("This protocol requires MiXCR 4.7.0.")
        root.mkdir(parents=True, exist_ok=False)
        failures = 0
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(execute_library, item, args, root, mixcr, trim): item[0]["analysis_library_id"] for item in items}
            for future in as_completed(futures):
                ident = futures[future]
                try:
                    future.result()
                    print(ident + ": completed", flush=True)
                except Exception as exc:
                    failures += 1
                    print(f"{ident}: {exc}", file=sys.stderr, flush=True)
        if failures:
            print(str(failures) + " libraries failed. Processing is incomplete.", file=sys.stderr)
            return 1
        print("MiXCR results written to the results directory.")
        return 0
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        print("ERROR: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
