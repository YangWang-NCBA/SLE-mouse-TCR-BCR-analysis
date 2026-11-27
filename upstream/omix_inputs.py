"""Build the analysis inputs from the OMIX downloads and the study's derived tables."""
from contextlib import ExitStack
import json
from pathlib import Path
import shutil
import zipfile

import pandas as pd

from fastqc_tables import extract_tables
from sequence_metrics import sample_id_from_library

METADATA_ZIP = "SLE_metadata_QC_PRJCA075761.zip"
FASTQC_ZIP = "SLE_FastQC_reports_PRJCA075761.zip"
LABELS = ["biological_sample_id", "mouse_id", "group_public", "tissue", "nominal_timepoint",
          "receptor", "analysis_library_id"]
OMITTED_INPUT_COLUMNS = {"study_day", "study_day_x", "study_day_y"}


def sources(omix, fastqc=None):
    """Accept the deposited ZIPs or extracted version 2.0 and legacy directories."""
    folder = omix.parent if omix.is_file() else omix
    metadata = omix if omix.is_file() else folder / METADATA_ZIP
    if not metadata.is_file():
        metadata = folder
        metadata_root(metadata)
    if fastqc is None:
        fastqc = folder / FASTQC_ZIP
        if not fastqc.is_file():
            fastqc = folder
    if not (fastqc.is_file() or (fastqc / "fastqc/raw").is_dir() or (fastqc / "raw").is_dir()):
        raise FileNotFoundError(f"Expected {FASTQC_ZIP} or the extracted FastQC reports. Use --fastqc to supply their location.")
    return metadata, fastqc


def read_table(path):
    compression = "gzip" if path.name.endswith(".gz") else None
    with path.open("rb") as handle:
        return pd.read_csv(handle, sep="\t", dtype=str, keep_default_na=False, compression=compression,
                           usecols=lambda name: name not in OMITTED_INPUT_COLUMNS)


def write_table(frame, path):
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def specimen_labels(libraries):
    pairs = libraries[["biological_sample_id", "analysis_library_id"]].copy()
    pairs["SampleID"] = pairs.analysis_library_id.map(sample_id_from_library)
    pairs = pairs[["biological_sample_id", "SampleID"]].drop_duplicates()
    return pairs.set_index("biological_sample_id", verify_integrity=True).SampleID


def with_sample_ids(frame, labels):
    if "analysis_library_id" in frame:
        library_ids = frame.analysis_library_id.unique()
        by_library = {key: sample_id_from_library(key) for key in library_ids}
        ids = frame.analysis_library_id.map(by_library)
    else:
        ids = frame.biological_sample_id.map(labels)
    if ids.isna().any():
        raise ValueError("A specimen has no matching sample label in the library metadata.")
    if "biological_sample_id" in frame and not ids.eq(frame.biological_sample_id.map(labels)).all():
        raise ValueError("The specimen and library identifiers refer to different sample labels.")
    for name in ("canonical_sample_label", "SampleID"):
        if name in frame and not frame[name].eq(ids).all():
            raise ValueError(f"{name} disagrees with the sample label in the library metadata.")
    result = frame.drop(columns=["biological_sample_id", "canonical_sample_label", "SampleID"], errors="ignore")
    result.insert(0, "SampleID", ids)
    return result


def prepared_input(source, scratch):
    """Use prepared inputs directly, adapting earlier table schemas when needed."""
    sample_file = source / "metadata/samples.tsv"
    if not sample_file.exists():
        sample_file = sample_file.with_suffix(".tsv.gz")
    headers = {path: set(pd.read_csv(path, sep="\t", nrows=0).columns)
               for path in source.rglob("*") if path.is_file() and path.name.endswith((".tsv", ".tsv.gz"))}
    migrate_ids = "SampleID" not in headers[sample_file]
    if not migrate_ids and all(OMITTED_INPUT_COLUMNS.isdisjoint(fields) for fields in headers.values()):
        return source
    if migrate_ids:
        library_file = source / "metadata/libraries.tsv"
        if not library_file.exists():
            library_file = library_file.with_suffix(".tsv.gz")
        labels = specimen_labels(read_table(library_file))
    for path in source.rglob("*"):
        if not path.is_file():
            continue
        target = scratch / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        fields = headers.get(path, set())
        needs_ids = migrate_ids and bool(fields.intersection({"analysis_library_id", "biological_sample_id"}))
        if needs_ids or fields.intersection(OMITTED_INPUT_COLUMNS):
            frame = read_table(path)
            if needs_ids:
                frame = with_sample_ids(frame, labels)
            write_table(frame, target)
        else:
            shutil.copyfile(path, target)
    return scratch


def join(left, right, keys, *, validate):
    merged = left.merge(right, on=keys, how="left", validate=validate, indicator=True)
    if not merged["_merge"].eq("both").all():
        raise ValueError(f"Missing OMIX or derived records for join on {keys}")
    return merged.drop(columns="_merge")


def metadata_root(root):
    """Locate a root layout or the one named wrapper used by earlier deposits."""
    stem = Path(METADATA_ZIP).stem
    for candidate in (root, root / stem, root / stem / stem):
        if (candidate / "samples.tsv").is_file() or (candidate / "metadata/sample_metadata.tsv").is_file():
            return candidate
    raise FileNotFoundError(f"Expected samples.tsv or metadata/sample_metadata.tsv under {root}")


def metadata_tables(root):
    """Restore the seven existing logical inputs from either deposited schema."""
    root = metadata_root(root)
    if not (root / "samples.tsv").is_file():
        meta, qc = root / "metadata", root / "qc_summary"
        tables = {name: read_table(meta / (name + ".tsv")) for name in (
            "sample_metadata", "library_metadata", "MiXCR_875_expected_chain_file_mapping",
            "native_TSV_851_manifest", "primer_sequences")}
        tables.update({name: read_table(qc / (name + ".tsv")) for name in (
            "mixcr_library_qc", "raw_fastq_library_qc")})
        if "experimental_week" not in tables["sample_metadata"]:
            crosswalk = read_table(meta / "sample_id_crosswalk.tsv")
            tables["sample_metadata"] = join(tables["sample_metadata"],
                crosswalk[["sample_id", "experimental_week"]], "sample_id", validate="one_to_one")
        return tables

    study_path = Path(__file__).resolve().parents[1] / "study.json"
    study = json.loads(study_path.read_text(encoding="utf-8"))
    if study.get("schema_version") != "2.0":
        raise ValueError(f"Unsupported metadata schema_version: {study.get('schema_version')!r}")
    samples = read_table(root / "samples.tsv").rename(columns={
        "week": "experimental_week", "time_label": "timepoint_label", "visit_order": "collection_order"})
    if not samples.timepoint_label.eq("W" + samples.experimental_week).all():
        raise ValueError("Sample collection weeks disagree with time labels")
    records = read_table(root / "libraries.tsv").sort_values(["sample_id", "library_id"])
    fastq = read_table(root / "fastq_qc.tsv")
    if not set(fastq.read_direction).issubset({"R1", "R2"}):
        raise ValueError("FASTQ read_direction must be R1 or R2")
    mates = fastq.pivot(index="library_id", columns="read_direction", values="filename")
    if list(mates.columns) != ["R1", "R2"] or mates.isna().any().any():
        raise ValueError("Each library must have both R1 and R2 FASTQ records")
    libraries = join(records[["library_id", "sample_id", "receptor"]],
        mates.rename(columns={"R1": "R1_public_filename", "R2": "R2_public_filename"}).reset_index(),
        "library_id", validate="one_to_one")
    annotations = join(records, samples[["sample_id", "mouse_id", "group", "tissue", "timepoint_label"]],
        "sample_id", validate="many_to_one").rename(columns={"sample_id": "biological_sample_id",
        "library_id": "analysis_library_id", "group": "group_public", "timepoint_label": "nominal_timepoint"})

    mixcr = annotations.copy()
    for target, source in (("tag_total", "align_total_reads_processed"),
                           ("refine_input_records", "align_aligned"), ("groups_total", "umi_output_diversity")):
        mixcr[target] = mixcr[source]
    mixcr["align_not_aligned"] = (mixcr.align_total_reads_processed.astype("int64")
                                  - mixcr.align_aligned.astype("int64")).astype(str)
    for target, numerator, denominator in (
        ("align_rate", "align_aligned", "align_total_reads_processed"),
        ("tag_match_rate", "tag_matched", "tag_total"),
        ("refine_output_rate", "refine_output_records", "refine_input_records"),
        ("assemble_reads_used_rate_of_refined", "assemble_reads_used", "refine_output_records")):
        mixcr[target] = [format(int(n) / int(d), ".10g") if int(d) else "NA"
                         for n, d in zip(mixcr[numerator], mixcr[denominator])]
    mixcr = mixcr[LABELS + ["align_total_reads_processed", "align_aligned", "align_rate", "align_not_aligned",
        "overlapped", "overlapped_aligned", "tag_total", "tag_matched", "tag_match_rate", "refine_input_records",
        "refine_output_records", "refine_output_rate", "umi_input_diversity", "umi_output_diversity",
        "umi_filtered_by_tag_quality", "groups_total", "groups_accepted", "assemble_final_clonotype_count",
        "assemble_reads_used", "assemble_aligned_reads_processed", "assemble_reads_used_rate_of_refined",
        "assemble_chain_count_closure_verified"]]

    counts = fastq.pivot(index="library_id", columns="read_direction", values="total_sequences").astype("int64")
    lengths = study["library_defaults"]
    raw = annotations.copy()
    raw["paired_read_count"] = raw.analysis_library_id.map(counts.R1).astype(str)
    raw["read_count_difference_R1_minus_R2"] = raw.analysis_library_id.map(counts.R1 - counts.R2).astype(str)
    bases = counts.R1 * int(lengths["read_length_R1_bp"]) + counts.R2 * int(lengths["read_length_R2_bp"])
    raw["total_bases_R1_plus_R2"] = raw.analysis_library_id.map(bases).map(lambda value: format(value, ".10g"))
    for name in ("file_qc_complete", "mean_read_length_bp"):
        raw[name] = study["observed_constants"]["raw_fastq_library_qc"][name]
    raw = raw[LABELS + ["file_qc_complete", "paired_read_count", "read_count_difference_R1_minus_R2",
        "total_bases_R1_plus_R2", "mean_read_length_bp", "weighted_Q20_percent", "weighted_Q30_percent",
        "weighted_mean_quality", "weighted_GC_percent"]]

    chain_records = read_table(root / "chain_records.tsv")
    if chain_records.duplicated(["library_id", "chain"]).any():
        raise ValueError("Repeated library/chain key in chain_records.tsv")
    states = {"present", "observed_zero_no_export", "not_targeted_export"}
    if not set(chain_records.file_state).issubset(states):
        raise ValueError("Unknown file_state in chain_records.tsv")
    chains = chain_records.loc[chain_records.file_state.ne("not_targeted_export")].copy()
    chains = join(chains, annotations[[*LABELS, "assemble_chain_count_closure_verified"]].rename(columns={
        "analysis_library_id": "library_id", "biological_sample_id": "sample_id", "group_public": "group"}),
        "library_id", validate="many_to_one").sort_values(["sample_id", "library_id", "chain"])
    chains["observed_file_count"] = chains.file_state.map({"present": "1", "observed_zero_no_export": "0"})
    chains["assemble_chain_functional"] = [str(int(total) - int(nonfunctional))
        if total not in ("", "NA") and nonfunctional not in ("", "NA") else "NA"
        for total, nonfunctional in zip(chains.assemble_chain_total, chains.assemble_chain_nonfunctional)]
    chains = chains.rename(columns={"expected_tsv_filename": "expected_public_filename",
        "assemble_chain_count_source": "count_source", "assemble_chain_count_closure_verified": "count_closure_verified"})
    status = chains[["sample_id", "mouse_id", "group", "tissue", "receptor", "library_id", "chain",
        "expected_public_filename", "file_state", "observed_file_count", "assemble_chain_total",
        "assemble_chain_nonfunctional", "assemble_chain_functional", "count_source", "count_closure_verified"]]
    native = chain_records.loc[chain_records.file_state.eq("present"),
        ["expected_tsv_filename", "library_id", "bytes", "data_rows"]].rename(columns={"expected_tsv_filename": "file_name"})
    return {"sample_metadata": samples, "library_metadata": libraries,
        "MiXCR_875_expected_chain_file_mapping": status, "native_TSV_851_manifest": native,
        "mixcr_library_qc": mixcr, "raw_fastq_library_qc": raw,
        "primer_sequences": read_table(root / "primers.tsv")}


def publication_tables(author, destination, samples, libraries, primers):
    tables = json.loads((author / "tables_1_to_3.json").read_text(encoding="utf-8"))
    sequenced = samples.loc[samples.biological_sample_id.isin(libraries.biological_sample_id)]
    for row in tables[0]["rows"][1:]:
        source, time = row[:2]
        for records in (samples, sequenced):
            subset = records if source == "Total" else records.loc[
                records.tissue.eq({"Peripheral blood": "blood", "Lymph nodes": "lymph_node"}.get(source, source.lower()))
                & records.nominal_week.eq(time.split()[-1])]
            counts = subset.groupby("group_public").biological_sample_id.nunique()
            row.append(" / ".join(str(counts.get(group, "—")) for group in ("IMQ_treated", "untreated_control")))

    primer_index = primers.set_index("primer_name", verify_integrity=True)
    used = []
    for row in tables[1]["rows"][1:]:
        for column, stage in ((2, "reverse_transcription"), (3, "PCR")):
            entries = []
            for name in row[column]:
                primer = primer_index.loc[name]
                if primer.stage != stage:
                    raise ValueError(f"Unexpected primer stage for {name}")
                used.append(name)
                entries += [name, primer.sequence_5_to_3]
            row[column] = "\n".join(entries)
    if sorted(used) != sorted(primer_index.index):
        raise ValueError("Table 2 primer names do not match the deposited primer table")
    (destination / "tables_1_to_3.json").write_text(json.dumps(tables, ensure_ascii=False, indent=2), encoding="utf-8")


def prepare_inputs(metadata, fastqc, author, destination):
    if fastqc is None:
        raise ValueError("Provide the FastQC archive or the directory containing the extracted FastQC reports.")
    for subdir in ("metadata", "analysis", "fastqc", "provenance"):
        (destination / subdir).mkdir(parents=True, exist_ok=True)
    with ExitStack() as opened:
        archive = zipfile.Path(opened.enter_context(zipfile.ZipFile(metadata))) if metadata.is_file() else metadata
        tables = metadata_tables(archive)
        samples = tables["sample_metadata"].rename(columns={
            "sample_id": "biological_sample_id", "group": "group_public", "collection_order": "visit_order",
            "timepoint_label": "nominal_timepoint", "experimental_week": "nominal_week",
        })
        samples = samples[["biological_sample_id", "mouse_id", "group_public", "tissue",
                           "nominal_timepoint", "nominal_week", "visit_order"]].copy()
        samples["group_code"] = samples.group_public.map({"untreated_control": "CON", "IMQ_treated": "IMQ"})
        if samples.group_code.isna().any() or not samples.biological_sample_id.is_unique:
            raise ValueError("Invalid group or repeated specimen ID in the OMIX sample table")
        libraries = tables["library_metadata"].rename(columns={
            "library_id": "analysis_library_id", "sample_id": "biological_sample_id",
            "R1_public_filename": "R1_filename", "R2_public_filename": "R2_filename",
        })[["biological_sample_id", "analysis_library_id", "receptor", "R1_filename", "R2_filename"]]
        libraries = join(libraries, samples, "biological_sample_id", validate="many_to_one")
        libraries = libraries.sort_values(["biological_sample_id", "analysis_library_id"])
        if not libraries.analysis_library_id.is_unique:
            raise ValueError("Repeated library ID in the OMIX library table")
        labels = specimen_labels(libraries)
        sample_table = with_sample_ids(samples, labels)
        if sample_table.duplicated(["SampleID", "tissue"]).any():
            raise ValueError("SampleID and tissue must identify one biological specimen.")
        write_table(sample_table, destination / "metadata/samples.tsv")
        public_libraries = with_sample_ids(libraries, labels)
        write_table(public_libraries, destination / "metadata/libraries.tsv")

        status = tables["MiXCR_875_expected_chain_file_mapping"].rename(columns={
            "sample_id": "biological_sample_id", "group": "group_public", "library_id": "analysis_library_id",
            "expected_public_filename": "expected_tsv_filename", "file_state": "chain_file_status",
            "count_source": "assemble_chain_count_source", "count_closure_verified": "assemble_chain_count_closure_verified",
        })
        status = join(status, samples[["biological_sample_id", "group_code", "nominal_timepoint"]],
                      "biological_sample_id", validate="many_to_one")
        write_table(with_sample_ids(status, labels), destination / "metadata/chains.tsv")
        # Keep file counts and byte sizes attached to the deposited native export they describe.
        native = tables["native_TSV_851_manifest"].rename(columns={
            "file_name": "expected_tsv_filename", "library_id": "analysis_library_id", "bytes": "file_size_bytes",
            "data_rows": "native_mixcr_rows_raw_file",
        })
        file_values = ["expected_tsv_filename", "analysis_library_id", "file_size_bytes", "native_mixcr_rows_raw_file"]
        present = status.loc[status.chain_file_status.eq("present"), ["analysis_library_id", "chain", "expected_tsv_filename"]]
        files = join(present, native[file_values], ["expected_tsv_filename", "analysis_library_id"], validate="one_to_one")
        files = files.rename(columns={"expected_tsv_filename": "tsv_filename"})
        analysis = destination / "analysis"
        metrics = read_table(author / "analysis/repertoire_chain_metrics.tsv.gz")
        metrics = join(metrics, files, ["analysis_library_id", "chain"], validate="one_to_one")
        metrics = join(metrics, libraries[LABELS], "analysis_library_id", validate="many_to_one")
        leading = LABELS + ["chain", "processing_status", "tsv_filename", "file_size_bytes", "native_mixcr_rows_raw_file"]
        metrics = metrics[leading + [c for c in metrics if c not in leading]]
        write_table(with_sample_ids(metrics, labels), analysis / "repertoire_chain_metrics.tsv")
        ranks = read_table(author / "analysis/rank_abundance_top100.tsv.gz")
        ranks = join(ranks, libraries[LABELS], "analysis_library_id", validate="many_to_one")
        ranks = with_sample_ids(ranks, labels)
        rank_columns = ["SampleID", *LABELS[1:], "chain", "rank", "clone_umi", "relative_umi_abundance"]
        write_table(ranks[rank_columns], analysis / "rank_abundance_top100.tsv")

        write_table(with_sample_ids(tables["mixcr_library_qc"], labels), analysis / "mixcr_library_qc.tsv")
        raw_qc = tables["raw_fastq_library_qc"]
        q30 = read_table(author / "analysis/read_q30.tsv")
        combined = join(q30[["analysis_library_id"]], raw_qc, "analysis_library_id", validate="one_to_one")
        combined = join(combined, q30, "analysis_library_id", validate="one_to_one")
        write_table(with_sample_ids(combined, labels), analysis / "raw_fastq_library_qc.tsv")
        publication_tables(author, destination, samples, libraries, tables["primer_sequences"])

    for name in ("analysis/native_counts.tsv", "provenance/unassigned_native_clonotypes.tsv.gz"):
        shutil.copyfile(author / name, destination / name)
    extract_tables(fastqc, public_libraries.to_dict("records"), destination / "fastqc")
