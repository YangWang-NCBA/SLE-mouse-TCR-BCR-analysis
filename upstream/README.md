# Generate the source tables

Run from the directory containing `reproduce.py`. The examples use `../sle-run/` for downloads and results. First run `python reproduce.py --omix ../sle-run/omix --output ../sle-run/results` with the two OMIX archives listed in the [main README](../README.md).

Install MiXCR 4.7.0, activate its licence and install the imgt.202312-3.sv8 reference in `~/.mixcr/libraries/`. Raw-read processing also uses Trim Galore, Cutadapt, FastQC (0.12.1 in this study) and SeqKit (2.14.0). The original Trim Galore and Cutadapt versions were not recorded. The `study.json` file at the root of the analysis package records the species, strain and supplier of the BALB/c mice, the SURFSeq 5000 instrument, the `N{4}TN{4}TN{4}` UMI pattern, Q20 trimming threshold, 35 bp minimum read length, and the full MiXCR build and tag pattern. Sample and library records are in `samples.tsv` and `libraries.tsv`; `dictionary.tsv` defines the consolidated fields. The input reader also supports earlier metadata packages and produces the same prepared tables used by the commands below. [Metadata package guide](../METADATA_PACKAGE.md) [MiXCR release](https://github.com/milaboratory/mixcr/releases/tag/v4.7.0) · [Licence activation](https://mixcr.com/mixcr/getting-started/milm/) · [Reference file](https://github.com/repseqio/library-imgt/releases/download/v8/imgt.202312-3.sv8.json.gz) · [Library guide](https://mixcr.com/mixcr/guides/external-libraries/)

### From raw FASTQ files

Place the paired reads listed in `../sle-run/results/data/metadata/libraries.tsv` under `../sle-run/fastq/`. Use a new directory for MiXCR output. Add `--library-id ID` to the first command to process one library.

```bash
python upstream/mixcr_processing.py --manifest ../sle-run/results/data/metadata/libraries.tsv --raw-root ../sle-run/fastq --out-root ../sle-run/processed --threads 8 --jobs 1
python upstream/sequence_metrics.py seqkit --metadata-dir ../sle-run/results/data/metadata --raw-root ../sle-run/fastq --out-root ../sle-run/metrics --threads 8
python upstream/sequence_metrics.py mixcr --metadata-dir ../sle-run/results/data/metadata --mixcr-root ../sle-run/processed/results --out-root ../sle-run/metrics
python upstream/sequence_metrics.py repertoire --metadata-dir ../sle-run/results/data/metadata --mixcr-root ../sle-run/processed/results --out-root ../sle-run/metrics --threads 8
```

The three metric commands produce read-quality, UMI, CDR3 and rank–abundance tables in `../sle-run/metrics/data/analysis/`. Run `mixcr` before `repertoire` so the latter can distinguish zero clonotypes from a missing export.

### From processed OMIX files

Extract the native MiXCR TSVs under `../sle-run/deposited_tsv/` to count distinct clone IDs within each library. The second command reads the deposited FastQC archive and produces the three input tables for Fig. 2.

```bash
python upstream/count_native_exports.py --tsv-root ../sle-run/deposited_tsv --release-dir ../sle-run/results --output-dir ../sle-run/native_counts
python upstream/fastqc_tables.py --fastqc-source ../sle-run/omix/SLE_FastQC_reports_PRJCA075761.zip --libraries ../sle-run/results/data/metadata/libraries.tsv --output ../sle-run/fastqc
```

### Use the regenerated tables

Copy `../sle-run/results/data/` to `../sle-run/recomputed-data/`. The table maps the example outputs to the files they replace in this copy.

| Generated file(s) | Destination in the copy |
|---|---|
| Three FastQC TSVs | `fastqc/` |
| `raw_fastq_library_qc.tsv`, `mixcr_library_qc.tsv`, `repertoire_chain_metrics.tsv`, `rank_abundance_top100.tsv` | `analysis/` with the same filenames |
| `mixcr_chain_status.tsv` | `metadata/chains.tsv` |
| `native_assembled_counts_by_library.tsv` | `analysis/native_counts.tsv` |
| `unassigned_native_clonotypes.tsv.gz` | `provenance/unassigned_native_clonotypes.tsv.gz` |

```bash
python reproduce.py --data-dir ../sle-run/recomputed-data --output ../sle-run/new-figures
```

## Data columns

`SampleID` is the label recorded at collection (e.g. `CON_8_1`). Use `SampleID` with `tissue` to identify one specimen. The same endpoint label appears in multiple organs from the same mouse. For control blood, labels such as `CON_1_1` and `CON_2_1` denote the week 0 and week 8 collections from mouse 1. Use `nominal_week` for collection week. In the phenotype table, `SampleID` matches the spleen photograph.

`mouse_id` identifies an animal across time and tissues. To retrieve the tissue for a packaged record, join its `analysis_library_id` to the same column in `results/data/metadata/libraries.tsv`. This library key matches `library_id` in the deposited metadata, which also lists the deposited specimen key (`sample_id`). `chain` is the receptor locus.

<details>
<summary>Measurements in the specimen–chain metrics table</summary>

Each row describes one native chain export from one specimen. An eligible record passes the chain and CDR3 checks, has a V and J call and a positive UMI count, and has no explicit nonproductive flag. Diversity, rank and length statistics use eligible records merged by V gene, J gene and nucleotide CDR3. Fractions range from 0 to 1. Missing or undefined measurements are empty or `NA`.

| Column | Definition |
|---|---|
| `SampleID`, `analysis_library_id`, `chain` | Collection label, sequencing library and receptor locus, as defined above. |
| `processing_status` | Processing outcome for the chain export (`ok` for a completed export). |
| `native_mixcr_rows_chain_valid` | Number of native export records assigned to the specified chain by the analysis rule. |
| `raw_file_total_umi` | Sum of UMIs across all records in this export. |
| `total_umi_chain_valid` | Sum of UMIs across chain-valid records. |
| `total_read_count` | Sum of available MiXCR read counts across all export records. |
| `analysis_eligible_rows` | Number of eligible records before merging. |
| `eligible_umi` | Sum of UMIs in eligible records. |
| `eligible_umi_fraction` | Eligible UMIs divided by all chain-valid UMIs. |
| `strict_clonotypes` | Number of distinct V gene, J gene and nucleotide CDR3 combinations among eligible records. Figs. 4a and 4c use native clonotype counts. |
| `shannon_entropy` | Shannon entropy, −Σ p ln(p), where p is the UMI fraction of a merged clonotype. |
| `shannon_effective` | Exponential of Shannon entropy, in effective clonotypes. |
| `simpson_effective` | Inverse Simpson index, 1/Σ p², in effective clonotypes. |
| `normalized_shannon_clonality` | One minus Shannon entropy divided by ln(number of merged clonotypes), defined for at least two clonotypes. |
| `top1_fraction`, `top10_fraction`, `top100_fraction` | Fraction of eligible UMIs carried by the first 1, 10 or 100 merged clonotypes ranked by UMI count. |
| `umi_weighted_mean_cdr3_aa_length` | Mean CDR3 length in amino acids, weighted by eligible UMI counts. |
| `median_minQualCDR3` | Median of the per-record minimum CDR3 base quality (Phred), over eligible records with a numeric MiXCR `minQualCDR3` value. |
| `missing_v_rows`, `missing_j_rows` | Export records without an assignable V or J gene, respectively. |
| `missing_cdr3_rows`, `cdr3_present_rows` | Counts of export records missing either CDR3 sequence, or containing both nucleotide and amino acid sequences, respectively. |
| `invalid_cdr3_rows` | Records with both CDR3 sequences that fail the translation or amino acid check. |
| `invalid_cdr3_fraction_of_present` | Invalid CDR3 records divided by records containing both sequences. |
| `v_top_tie_rows`, `j_top_tie_rows` | Export records with distinct V or J genes tied for the highest alignment score. |
| `chain_conflict_rows` | Export records that fail the assignment rule for this chain. |
| `explicit_productive_false_rows` | Chain-valid records explicitly marked nonproductive in the native export. |
| `error_message` | Processing error message, empty or `NA` for successful records. |

</details>

<details>
<summary>Native counts, rank–abundance, read quality and figure summaries</summary>

| Table | Columns |
|---|---|
| `native_counts.tsv` | `export_files` counts chain exports per library. `distinct_assembled_clone_ids` counts distinct native clone IDs across them. `duplicate_rows_across_exports` counts repeat appearances of the same ID. `IGH`, `IGK`, `IGL`, `TRA`, `TRB`, `TRD` and `TRG` count IDs assigned to each locus. `unassigned_chain` counts the remaining IDs. |
| `unassigned_native_clonotypes.tsv.gz` | `cloneId` is the native MiXCR clone ID within its library. `chains` is the original chain assignment, including empty entries. |
| `rank_abundance_top100.tsv.gz` | `rank` starts at 1 for the most abundant merged clonotype within a specimen and chain. `clone_umi` is its UMI count. `relative_umi_abundance` is its fraction of all eligible UMIs, including those ranked below 100. |
| `read_q30.tsv` | `R1_Q30_percent` and `R2_Q30_percent` give the percentage of bases with Phred quality ≥30 in each raw mate file. These values support the median Q30 base percentages and ranges reported separately for R1 and R2 in the manuscript. |
| Fig. 2 and Fig. 4d summaries | `n` counts contributing files at a base position (Fig. 2a,c) or contributing specimen–chain measurements at a rank (Fig. 4d). `q25` and `q75` are the 25th and 75th percentiles. |
| Fig. 3 summary | `lower`, `middle` and `upper` are the 25th, 50th and 75th percentiles of library output UMI counts. `ymin` and `ymax` are the lowest and highest observations within 1.5 interquartile ranges of the box. |

</details>
