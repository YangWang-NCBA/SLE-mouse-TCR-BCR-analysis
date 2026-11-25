"""Parse MiXCR reports and calculate the manuscript sequence metrics."""


from __future__ import annotations
import argparse, concurrent.futures as cf, csv, json, math, re, statistics, sys, subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any,Iterable


NA = "NA"


def scan_files(root, predicate):
    return sorted((p for p in root.rglob("*") if p.is_file() and predicate(p)), key=lambda p: (p.name,str(p)))


CHAINS = {"BCR": ("IGH", "IGK", "IGL"), "TCR": ("TRA", "TRB", "TRD", "TRG")}


CHAIN_SET = set(sum((list(x) for x in CHAINS.values()), []))


CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def sample_id_from_library(library: str) -> str:
    match = re.fullmatch(r".+_((?:CON|IMQ)_\d+_\d+)_(?:BCR|TCR)_L00", library)
    if not match:
        raise ValueError(f"Unrecognized study library ID: {library}")
    return match[1]


def clean_value(value: Any) -> str:
    if value is None:
        return NA
    if isinstance(value, float):
        if not math.isfinite(value):
            return NA
        return format(value, ".10g")
    value = str(value)
    return NA if value == "" else value.replace("\r", " ").replace("\n", " ")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and "analysis_library_id" in rows[0] and "SampleID" not in rows[0]:
        for row in rows:
            row["SampleID"] = sample_id_from_library(row["analysis_library_id"])
    return rows


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str] | None = None) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_value(row.get(key)) for key in fields})


def fnum(value: Any) -> float | None:
    if value in (None, "", NA):
        return None
    try:
        result = float(str(value).replace(",", ""))
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def fint(value: Any) -> int | None:
    number = fnum(value)
    return None if number is None else int(round(number))


def ratio(num: Any, den: Any) -> float | None:
    a, b = fnum(num), fnum(den)
    return None if a is None or b in (None, 0) else a / b


def metadata_index(metadata_dir: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    libraries = read_tsv(metadata_dir / "libraries.tsv")
    by_library = {r["analysis_library_id"]: r for r in libraries}
    if len(by_library) != len(libraries):
        raise ValueError("analysis_library_id is not unique in libraries.tsv")
    return libraries, by_library


def parse_fastqc_modules(text: str) -> list[tuple[str, str, list[str]]]:
    modules, name, status, lines = [], None, None, []
    for line in text.splitlines():
        if line.startswith(">>") and line != ">>END_MODULE":
            parts = line[2:].split("\t")
            name, status, lines = parts[0], parts[1] if len(parts) > 1 else NA, []
        elif line == ">>END_MODULE":
            if name is not None:
                modules.append((name, status or NA, lines))
            name, status, lines = None, None, []
        elif name is not None:
            lines.append(line)
    return modules


def parse_table(lines: list[str]) -> tuple[list[str], list[list[str]]]:
    header, rows = [], []
    for line in lines:
        if line.startswith("#") and not header:
            header = line[1:].split("\t")
        elif line and not line.startswith("#"):
            rows.append(line.split("\t"))
    return header, rows


def library_from_report_name(name: str) -> tuple[str | None, str | None]:
    patterns = [
        (r"^(.*)_align\.report\.json$", "align_json"), (r"^(.*)_align\.report\.txt$", "align_text"),
        (r"^(.*)\.refine\.report\.txt$", "refine_text"), (r"^(.*)\.assemble\.report\.txt$", "assemble_text"),
        (r"^(.*)\.findAlleles\.report\.txt$", "findAlleles_text"),
    ]
    for pattern, kind in patterns:
        match = re.match(pattern, name)
        if match:
            return match.group(1), kind
    return None, None


def text_metric(text: str, label: str) -> float | None:
    match = re.search(rf"^{re.escape(label)}:\s*([0-9][0-9,]*(?:\.[0-9]+)?)", text, re.M)
    return fnum(match.group(1)) if match else None


def assemble_chain_counts(text: str) -> tuple[dict[str, tuple[int | None, int | None, str]], bool]:
    """Infer omitted zero chains only from a completed, numerically closed report."""
    total_matches = list(re.finditer(r"^\s*([A-Z]+) chains:\s*([0-9][0-9,]*)", text, re.M))
    nf_matches = list(re.finditer(r"^\s*([A-Z]+) non-functional:\s*([0-9][0-9,]*)", text, re.M))
    total_labels = re.findall(r"^\s*([A-Z]+) chains:", text, re.M)
    nf_labels = re.findall(r"^\s*([A-Z]+) non-functional:", text, re.M)
    if len(total_matches) != len(total_labels) or len(nf_matches) != len(nf_labels):
        raise ValueError("invalid numeric chain counter in assemble report")
    totals = {m.group(1): fint(m.group(2)) for m in total_matches}
    nonfunctional = {m.group(1): fint(m.group(2)) for m in nf_matches}
    if len(totals) != len(total_matches) or len(nonfunctional) != len(nf_matches):
        raise ValueError("duplicate chain counter labels in assemble report")
    final_count = text_metric(text, "Final clonotype count")
    last_line = text.rstrip().splitlines()[-1] if text.rstrip() else ""
    footer_complete = bool(re.fullmatch(r"={10,}", last_line.strip()))
    count_closed = final_count is not None and final_count >= 0 and final_count == int(final_count) and sum(totals.values()) == final_count
    # A composite label such as TRAD leaves the individual chains unresolved.
    closure_verified = footer_complete and count_closed and set(totals).issubset(CHAIN_SET)
    result = {}
    for chain in CHAIN_SET:
        total = totals.get(chain, 0 if closure_verified else None)
        nf = nonfunctional.get(chain, 0 if total == 0 else None)
        if total is not None and nf is not None and nf > total:
            raise ValueError(f"assemble non-functional count exceeds total for {chain}")
        source = "explicit" if chain in totals else ("inferred_zero_closed_report" if closure_verified else "missing")
        result[chain] = (total, nf, source)
    return result, closure_verified


def cmd_mixcr(args: argparse.Namespace) -> None:
    mixcr_root, out_root, metadata_dir = Path(args.mixcr_root).resolve(), Path(args.out_root).resolve(), Path(args.metadata_dir).resolve()
    libraries, by_library = metadata_index(metadata_dir)
    report_files = scan_files(mixcr_root, lambda p: library_from_report_name(p.name)[0] is not None)
    by_key: dict[tuple[str, str], list[Path]] = defaultdict(list)
    for path in report_files:
        lib, kind = library_from_report_name(path.name)
        if lib in by_library and kind:
            by_key[(lib, kind)].append(path)
    qc_rows, chain_rows = [], []
    expected_kinds = ("align_json", "refine_text", "assemble_text")
    for lib in libraries:
        lib_id = lib["analysis_library_id"]
        for kind in expected_kinds:
            if len(by_key[(lib_id, kind)]) != 1:
                raise ValueError(f"Expected one {kind} report for {lib_id}")
        base = {k: lib[k] for k in ("SampleID", "mouse_id", "group_public", "tissue", "nominal_timepoint", "receptor", "analysis_library_id")}
        row: dict[str, Any] = dict(base)
        data = json.loads(by_key[(lib_id, "align_json")][0].read_text(encoding="utf-8"))
        row.update({
            "align_total_reads_processed": data.get("totalReadsProcessed"), "align_aligned": data.get("aligned"),
            "align_rate": ratio(data.get("aligned"), data.get("totalReadsProcessed")), "align_not_aligned": data.get("notAligned"),
            "overlapped": data.get("overlapped"), "overlapped_aligned": data.get("overlappedAligned"),
            "tag_total": (data.get("tagParsingReport") or {}).get("total"), "tag_matched": (data.get("tagParsingReport") or {}).get("matched"),
            "tag_match_rate": ratio((data.get("tagParsingReport") or {}).get("matched"), (data.get("tagParsingReport") or {}).get("total")),
        })
        usage = data.get("chainUsage") or {}
        for chain, values in (usage.get("chains") or {}).items():
            chain_rows.append({"SampleID": lib["SampleID"], "analysis_library_id": lib_id, "chain": chain, "align_chain_total": values.get("total"), "align_chain_nonfunctional": values.get("nonFunctional"), "align_chain_out_of_frame": values.get("isOOF"), "align_chain_has_stops": values.get("hasStops")})
        text = by_key[(lib_id, "refine_text")][0].read_text(encoding="utf-8", errors="replace")
        row.update({
            "refine_input_records": text_metric(text, "Number of input records"), "refine_output_records": text_metric(text, "Number of output records"),
            "refine_output_rate": ratio(text_metric(text, "Number of output records"), text_metric(text, "Number of input records")),
            "umi_input_diversity": text_metric(text, "  UMI input diversity"), "umi_output_diversity": text_metric(text, "  UMI output diversity"),
            "umi_filtered_by_tag_quality": text_metric(text, "  UMI diversity filtered by tag quality"),
            "groups_total": text_metric(text, "  Number of groups"), "groups_accepted": text_metric(text, "  Number of groups accepted"),
        })
        text = by_key[(lib_id, "assemble_text")][0].read_text(encoding="utf-8", errors="replace")
        row.update({
            "assemble_final_clonotype_count": text_metric(text, "Final clonotype count"), "assemble_reads_used": text_metric(text, "Reads used in clonotypes, percent of total"),
            "assemble_aligned_reads_processed": text_metric(text, "Aligned reads processed"),
            "assemble_reads_used_rate_of_refined": ratio(text_metric(text, "Reads used in clonotypes, percent of total"), row.get("refine_output_records")),
        })
        parsed_chains, closure_verified = assemble_chain_counts(text)
        row["assemble_chain_count_closure_verified"] = "yes" if closure_verified else "no"
        for chain, (total, nf, source) in parsed_chains.items():
            functional = None if total is None or nf is None else total - nf
            chain_rows.append({"SampleID": lib["SampleID"], "analysis_library_id": lib_id, "chain": chain, "assemble_chain_total": total, "assemble_chain_nonfunctional": nf, "assemble_chain_functional": functional, "assemble_chain_count_source": source, "assemble_chain_count_closure_verified": "yes" if closure_verified else "no", "assemble_report_observed": "yes"})
        required_metrics = ("align_total_reads_processed", "align_aligned", "refine_input_records", "refine_output_records", "assemble_final_clonotype_count")
        unparsed = [field for field in required_metrics if fnum(row.get(field)) is None]
        if unparsed:
            raise ValueError(f"Missing MiXCR metrics in {lib_id}: {', '.join(unparsed)}")
        qc_rows.append(row)
    write_tsv(out_root / "data/analysis/mixcr_library_qc.tsv", qc_rows)
    write_tsv(out_root / "data/analysis/mixcr_chain_report_counts_long.tsv", chain_rows)
def parse_hit(hit_string: str | None) -> tuple[str | None, bool]:
    if not hit_string or hit_string == NA:
        return None, False
    parsed = []
    for part in str(hit_string).split(","):
        match = re.match(r"^(.+)\((-?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\)$", part.strip())
        if match:
            parsed.append((match.group(1).strip(), float(match.group(2))))
    if not parsed:
        return None, False
    best_score = max(score for _, score in parsed)
    calls = [call for call, score in parsed if score == best_score]
    genes = []
    for call in calls:
        gene = re.sub(r"\*.*$", "", call)
        if gene not in genes:
            genes.append(gene)
    return genes[0] if genes else None, len(genes) > 1


def translate_nt(nt: str) -> str | None:
    nt = nt.strip().upper()
    if not nt or re.search(r"[^ACGT]", nt) or len(nt) % 3:
        return None
    aa = "".join(CODON.get(nt[i:i+3], "?") for i in range(0, len(nt), 3))
    return None if "?" in aa else aa


def chain_matches(raw: str, chain: str) -> bool:
    return bool(re.search(rf"(^|[,;|\s]){re.escape(chain)}($|[,;|\s])", (raw or "").upper()))


def _process_chain(task: tuple[dict[str, str], str]) -> tuple[dict[str, Any], list[dict[str, Any]], str | None]:
    meta, path_str = task
    path, chain = Path(path_str), meta["chain"]
    base = {k: meta.get(k, NA) for k in ("SampleID", "mouse_id", "group_public", "tissue", "nominal_timepoint", "receptor", "analysis_library_id", "chain")}
    required = {"uniqueMoleculeCount", "allVHitsWithScore", "allJHitsWithScore", "nSeqCDR3", "aaSeqCDR3", "chains"}
    native_raw = chain_valid_n = total_file_umi = total_umi = eligible_rows = eligible_umi = total_reads = 0.0
    missing_v = missing_j = missing_cdr3 = invalid_cdr3 = cdr3_present = chain_conflicts = productive_false_n = 0
    v_ties = j_ties = 0
    strict: dict[tuple[str, str, str], list[Any]] = {}
    minquals: list[float] = []
    try:
        csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = set(reader.fieldnames or [])
            missing = sorted(required - header)
            if missing:
                raise ValueError("missing required columns: " + ",".join(missing))
            productive_col = next((x for x in ("isProductive", "productive", "Productive") if x in header), None)
            for row in reader:
                native_raw += 1
                umi = fnum(row.get("uniqueMoleculeCount"))
                if umi is None or umi < 0 or abs(umi - round(umi)) > 1e-8:
                    raise ValueError(f"invalid uniqueMoleculeCount at data row {int(native_raw)}")
                total_file_umi += umi
                read_count = fnum(row.get("readCount"))
                if read_count is not None:
                    total_reads += read_count
                v_gene, v_tie = parse_hit(row.get("allVHitsWithScore"))
                j_gene, j_tie = parse_hit(row.get("allJHitsWithScore"))
                v_ties += int(v_tie); j_ties += int(j_tie)
                nt = (row.get("nSeqCDR3") or "").strip().upper()
                aa = (row.get("aaSeqCDR3") or "").strip().upper()
                present = bool(nt and aa)
                cdr3_present += int(present)
                translated = translate_nt(nt)
                valid_cdr3 = bool(present and translated == aa and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", aa))
                invalid_cdr3 += int(present and not valid_cdr3)
                if not v_gene: missing_v += 1
                if not j_gene: missing_j += 1
                if not present: missing_cdr3 += 1
                valid_chain = bool(j_gene and j_gene.startswith("TRDJ")) if chain == "TRD" else chain_matches(row.get("chains") or "", chain)
                if not valid_chain:
                    chain_conflicts += 1
                    continue
                chain_valid_n += 1
                total_umi += umi
                prod_false = False
                if productive_col:
                    prod_false = (row.get(productive_col) or "").strip().upper() in {"FALSE", "F", "0", "NO", "N"}
                productive_false_n += int(prod_false)
                eligible = umi > 0 and bool(v_gene) and bool(j_gene) and valid_cdr3 and not prod_false
                if eligible:
                    eligible_rows += 1
                    eligible_umi += umi
                    key = (v_gene or "", j_gene or "", nt)
                    if key in strict:
                        strict[key][0] += umi
                    else:
                        strict[key] = [umi, aa]
                    q = fnum(row.get("minQualCDR3"))
                    if q is not None:
                        minquals.append(q)
        weights = sorted((float(v[0]), key, str(v[1])) for key, v in strict.items())
        weights.reverse()
        total_weight = sum(w for w, _, _ in weights)
        probs = [w / total_weight for w, _, _ in weights] if total_weight > 0 else []
        richness = len(weights)
        entropy = -sum(p * math.log(p) for p in probs if p > 0)
        shannon_eff = math.exp(entropy) if probs else None
        simpson_eff = 1 / sum(p * p for p in probs) if probs else None
        clonality = 1 - entropy / math.log(richness) if richness > 1 else None
        mean_len = sum(len(aa) * w for w, _, aa in weights) / total_weight if total_weight else None
        def top_fraction(n: int) -> float | None:
            return sum(w for w, _, _ in weights[:n]) / total_weight if total_weight else None
        metric = {**base, "processing_status": "ok", "tsv_filename": path.name, "file_size_bytes": path.stat().st_size,
                  "native_mixcr_rows_raw_file": int(native_raw), "native_mixcr_rows_chain_valid": int(chain_valid_n),
                  "raw_file_total_umi": total_file_umi, "total_umi_chain_valid": total_umi, "total_read_count": total_reads,
                  "analysis_eligible_rows": int(eligible_rows), "eligible_umi": eligible_umi, "eligible_umi_fraction": ratio(eligible_umi, total_umi),
                  "strict_clonotypes": richness, "shannon_entropy": entropy if probs else None, "shannon_effective": shannon_eff,
                  "simpson_effective": simpson_eff, "normalized_shannon_clonality": clonality,
                  "top1_fraction": top_fraction(1), "top10_fraction": top_fraction(10), "top100_fraction": top_fraction(100),
                  "umi_weighted_mean_cdr3_aa_length": mean_len, "median_minQualCDR3": statistics.median(minquals) if minquals else None,
                  "missing_v_rows": missing_v, "missing_j_rows": missing_j, "missing_cdr3_rows": missing_cdr3,
                  "cdr3_present_rows": cdr3_present, "invalid_cdr3_rows": invalid_cdr3,
                  "invalid_cdr3_fraction_of_present": ratio(invalid_cdr3, cdr3_present), "v_top_tie_rows": v_ties, "j_top_tie_rows": j_ties,
                  "chain_conflict_rows": chain_conflicts, "explicit_productive_false_rows": productive_false_n, "error_message": NA}
        ranks = [{**base, "rank": rank, "clone_umi": w, "relative_umi_abundance": w / total_weight} for rank, (w, _, _) in enumerate(weights[:100], 1)]
        return metric, ranks, None
    except Exception as exc:
        error = str(exc).replace(path_str, path.name)
        return {**base, "processing_status": "error", "tsv_filename": path.name, "error_message": error}, [], error


def cmd_repertoire(args: argparse.Namespace) -> None:
    mixcr_root, out_root, metadata_dir = Path(args.mixcr_root).resolve(), Path(args.out_root).resolve(), Path(args.metadata_dir).resolve()
    expected = [{k: r[k] for k in ("SampleID", "mouse_id", "group_code", "group_public", "tissue", "nominal_timepoint", "receptor", "analysis_library_id", "chain", "expected_tsv_filename")} for r in read_tsv(metadata_dir / "chains.tsv")]
    chain_files = scan_files(mixcr_root, lambda p: bool(re.search(r"_(?:IGH|IGK|IGL|TRA|TRB|TRD|TRG)\.tsv$", p.name)))
    by_name: dict[str, list[Path]] = defaultdict(list)
    for path in chain_files:
        by_name[path.name].append(path)
    report_long = read_tsv(out_root / "data/analysis/mixcr_chain_report_counts_long.tsv")
    assemble: dict[tuple[str, str], dict[str, str]] = {}
    for row in report_long:
        if row.get("assemble_report_observed") == "yes":
            assemble[(row["analysis_library_id"], row["chain"])] = row
    status_rows, tasks = [], []
    for row in expected:
        paths = by_name.get(row["expected_tsv_filename"], [])
        rep = assemble.get((row["analysis_library_id"], row["chain"]))
        if len(paths) == 1:
            status = "present"
            tasks.append((row, str(paths[0])))
        elif len(paths) > 1:
            status = "duplicate_chain_file"
        elif rep is None:
            status = "report_missing"
        else:
            total = fint(rep.get("assemble_chain_total"))
            nonfunctional = fint(rep.get("assemble_chain_nonfunctional"))
            if total is None or nonfunctional is None:
                status = "report_metrics_missing"
            elif nonfunctional > total or total < 0 or nonfunctional < 0:
                status = "report_metrics_invalid"
            elif total == 0:
                closed = rep.get("assemble_chain_count_closure_verified") == "yes"
                status = "observed_zero_no_export" if closed else "report_metrics_invalid"
            else:
                status = "nonfunctional_only_no_export" if nonfunctional == total else "export_missing_review"
        status_rows.append({**row, "observed_file_count": len(paths), "chain_file_status": status,
                            "assemble_chain_total": rep.get("assemble_chain_total") if rep else None,
                            "assemble_chain_nonfunctional": rep.get("assemble_chain_nonfunctional") if rep else None,
                            "assemble_chain_functional": rep.get("assemble_chain_functional") if rep else None,
                            "assemble_chain_count_source": rep.get("assemble_chain_count_source") if rep else None,
                            "assemble_chain_count_closure_verified": rep.get("assemble_chain_count_closure_verified") if rep else None})
    write_tsv(out_root / "data/analysis/mixcr_chain_status.tsv", status_rows)
    duplicates = [r for r in status_rows if r["chain_file_status"] == "duplicate_chain_file"]
    if duplicates:
        raise RuntimeError(f"Duplicate expected chain TSV basenames detected for {len(duplicates)} slots; no duplicate files were summed")
    unresolved = [r for r in status_rows if r['chain_file_status'] not in ('present','observed_zero_no_export','nonfunctional_only_no_export')]
    if unresolved:
        raise RuntimeError(f"Unresolved chain exports: {len(unresolved)}; inspect mixcr_chain_status.tsv")
    metrics, ranks = [], []
    workers = max(1, min(args.threads, 8))
    with cf.ProcessPoolExecutor(max_workers=workers) as pool:
        for task, result in zip(tasks, pool.map(_process_chain, tasks)):
            metric, rank_rows, error = result
            if error:
                raise RuntimeError(f"{Path(task[1]).name}: {error}")
            metrics.append(metric)
            ranks.extend(rank_rows)
    metrics.sort(key=lambda r: (r.get("mouse_id", ""), r.get("tissue", ""), r.get("nominal_timepoint", ""), r.get("receptor", ""), r.get("chain", "")))
    ranks.sort(key=lambda r: (r.get("mouse_id", ""), r.get("tissue", ""), r.get("nominal_timepoint", ""), r.get("chain", ""), int(r.get("rank", 0))))
    write_tsv(out_root / "data/analysis/repertoire_chain_metrics.tsv", metrics)
    write_tsv(out_root / "data/analysis/rank_abundance_top100.tsv", ranks)
def norm_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())

def _seqkit_one(task: tuple[dict[str, str], str]) -> dict[str, Any]:
    meta, path = task
    cmd = ["seqkit", "stats", "-a", "-T", "-j", "1", path]
    proc = subprocess.run(cmd, text=True, capture_output=True)
    base = {k: meta.get(k, NA) for k in ("filename", "SampleID", "analysis_library_id", "read_direction")}
    if proc.returncode != 0:
        error = proc.stderr.strip().replace(path, meta.get("filename", Path(path).name))[-800:]
        return {**base, "status": "error", "error_message": error}
    rows = list(csv.DictReader(proc.stdout.splitlines(), delimiter="\t"))
    if not rows:
        return {**base, "status": "error", "error_message": "seqkit returned no data row"}
    raw = {norm_header(k): v for k, v in rows[-1].items()}
    aliases = {
        "num_reads": ("numseqs", "numsequences"), "total_bases": ("sumlen",),
        "min_length": ("minlen",), "mean_length": ("avglen", "meanlen"), "max_length": ("maxlen",),
        "q20_percent": ("q20",), "q30_percent": ("q30",), "mean_quality": ("avgqual", "meanqual"),
        "gc_percent": ("gc",), "n_bases": ("sumn",),
    }
    out = {**base, "status": "ok", "error_message": NA}
    for target, choices in aliases.items():
        out[target] = next((raw[c] for c in choices if c in raw), NA)
    return out

def weighted_mean(rows: list[dict[str, Any]], field: str, weight: str) -> float | None:
    pairs = [(fnum(r.get(field)), fnum(r.get(weight))) for r in rows]
    pairs = [(x, w) for x, w in pairs if x is not None and w is not None and w > 0]
    return None if not pairs else sum(x * w for x, w in pairs) / sum(w for _, w in pairs)

def cmd_seqkit(args: argparse.Namespace) -> None:
    out_root = Path(args.out_root).resolve()
    from mixcr_processing import load_inputs
    meta=Path(args.metadata_dir)
    inputs=load_inputs(meta/'libraries.tsv',Path(args.raw_root))
    tasks=[]
    for row,mates in inputs:
        if "SampleID" not in row:
            row["SampleID"] = sample_id_from_library(row["analysis_library_id"])
        for direction,path in zip(('R1','R2'),mates):
            tasks.append(({**row,'filename':path.name,'read_direction':direction},str(path)))
    with cf.ThreadPoolExecutor(max_workers=max(1, min(args.threads, 24))) as pool:
        results = list(pool.map(_seqkit_one, tasks))
    for row in results:
        if row['status'] != 'ok':
            raise RuntimeError(f"SeqKit failed for {row['filename']}: {row['error_message']}")

    libraries = [row for row,_ in inputs]
    by_lib: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in results:
        by_lib[row["analysis_library_id"]].append(row)
    aggregate = []
    for lib in libraries:
        rr = by_lib.get(lib["analysis_library_id"], [])
        ok = [r for r in rr if r.get("status") == "ok"]
        by_direction = {r.get("read_direction"): r for r in ok}
        r1_count = fnum(by_direction.get("R1", {}).get("num_reads"))
        r2_count = fnum(by_direction.get("R2", {}).get("num_reads"))
        read_counts = [x for x in (r1_count, r2_count) if x is not None]
        total_bases = sum(x for x in (fnum(r.get("total_bases")) for r in ok) if x is not None)
        aggregate.append({
            **{k: lib[k] for k in ("SampleID", "mouse_id", "group_public", "tissue", "nominal_timepoint", "receptor", "analysis_library_id")},
            "file_qc_complete": "yes" if len(ok) == 2 else "no",
            "paired_read_count": min(read_counts) if len(read_counts) == 2 else None,
            "read_count_difference_R1_minus_R2": (r1_count - r2_count) if r1_count is not None and r2_count is not None else None,
            "total_bases_R1_plus_R2": total_bases if ok else None,
            "mean_read_length_bp": weighted_mean(ok, "mean_length", "num_reads"),
            "weighted_Q20_percent": weighted_mean(ok, "q20_percent", "total_bases"),
            "weighted_Q30_percent": weighted_mean(ok, "q30_percent", "total_bases"),
            "R1_Q30_percent": by_direction.get("R1", {}).get("q30_percent"),
            "R2_Q30_percent": by_direction.get("R2", {}).get("q30_percent"),
            "weighted_mean_quality": weighted_mean(ok, "mean_quality", "total_bases"),
            "weighted_GC_percent": weighted_mean(ok, "gc_percent", "total_bases"),
        })
    write_tsv(out_root / "data/analysis/raw_fastq_library_qc.tsv", aggregate)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['mixcr','repertoire','seqkit'])
    parser.add_argument('--mixcr-root')
    parser.add_argument('--raw-root')
    parser.add_argument('--out-root',required=True)
    parser.add_argument('--metadata-dir',required=True,help='Prepared data/metadata directory from reproduce.py')
    parser.add_argument('--threads',type=int,default=8)
    args=parser.parse_args()
    if args.command=='seqkit' and not args.raw_root:parser.error('--raw-root is required for seqkit')
    if args.command!='seqkit' and not args.mixcr_root:parser.error('--mixcr-root is required')
    {'mixcr':cmd_mixcr,'repertoire':cmd_repertoire,'seqkit':cmd_seqkit}[args.command](args)

if __name__=='__main__':main()
