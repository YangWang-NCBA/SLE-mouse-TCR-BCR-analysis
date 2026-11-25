"""Reproduce Figures 2–4 and Tables 1–4."""

import argparse
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parent
CODE = ROOT / "code"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New results directory")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--omix", type=Path, help="OMIX downloads directory, metadata ZIP, or extracted directory")
    inputs.add_argument("--data-dir", type=Path, help="Prepared data directory from a previous run or upstream recalculation")
    parser.add_argument("--fastqc", type=Path, help="FastQC report ZIP or extracted directory if stored separately")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--r-library", type=Path, help="Additional R library directory")
    parser.add_argument("--formats", default="pdf,png", help="Include pdf,png and optionally svg,tiff")
    parser.add_argument("--font", choices=["Arial", "Liberation Sans"], default="Arial")
    args = parser.parse_args()

    formats = args.formats.split(",")
    selected = set(formats)
    if len(selected) != len(formats) or not {"pdf", "png"} <= selected <= {"pdf", "png", "svg", "tiff"}:
        parser.error("Include pdf,png once each. Optional formats are svg and tiff.")

    output = args.output.resolve()
    source = (args.data_dir or args.omix).resolve()
    if (output == ROOT or output in ROOT.parents or output == source or output in source.parents
            or ROOT / "data" in output.parents or (args.data_dir and source in output.parents)):
        parser.error("Choose a new directory separate from the input data and the package root.")
    if output.exists():
        parser.error("Choose a new output directory.")
    if args.fastqc and not args.omix:
        parser.error("--fastqc requires --omix.")
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT / "upstream"))
    from omix_inputs import prepare_inputs, prepared_input, sources
    if args.omix:
        try:
            metadata, fastqc = sources(source, args.fastqc.resolve() if args.fastqc else None)
        except FileNotFoundError as error:
            parser.error(str(error))
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as tmp:
        working = Path(tmp) / "results"
        working.mkdir()
        data = source
        if args.omix:
            data = working / "data"
            prepare_inputs(metadata, fastqc, ROOT / "data", data)
        else:
            data = prepared_input(source, Path(tmp) / "prepared")
        subprocess.run([sys.executable, "-B", CODE / "tables.py", data, working], check=True)
        plot = [args.rscript, "--vanilla", CODE / "plot.R", ROOT, working, args.formats, data]
        if args.r_library:
            plot.append(args.r_library.resolve())
        env = dict(os.environ, SD_FONT_FAMILY=args.font, SD_BOLD_FONT_FILE=str(Path(tmp) / "bold.txt"))
        subprocess.run(plot, check=True, env=env)
        subprocess.run([sys.executable, "-B", CODE / "assemble.py", working, args.formats], check=True, env=env)
        working.rename(output)

    print(f"Figures and tables written to {output.name}")


if __name__ == "__main__":
    main()
