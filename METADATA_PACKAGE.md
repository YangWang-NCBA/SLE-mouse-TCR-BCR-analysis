# SLE metadata and quality summaries

This guide describes the six data tables in `SLE_metadata_QC_PRJCA075761.zip` (OMIX022213). The analysis package includes `study.json` at its root. It records project accessions, shared experimental conditions, processing parameters, software versions, native file schemas and provenance. The metadata archive includes `dictionary.tsv` with field definitions, units, original aliases and recovery formulas. [Study parameters](study.json)

| OMIX metadata file | Records | Purpose |
|---|---:|---|
| samples.tsv | 125 | Biological specimens, mouse/group/tissue and sampling time |
| libraries.tsv | 250 | Specimen-to-library links and independent library QC values |
| fastq_qc.tsv | 500 | Original FASTQ filenames, MD5 values and raw-read FastQC summaries |
| chain_records.tsv | 875 | Expected chain-export slots and their QC/count information |
| primers.tsv | 23 | Reverse-transcription and PCR primers |
| dictionary.tsv | — | Field definitions and recovery rules |

Join samples to libraries by `sample_id`, and libraries to FASTQ/chain records by `library_id`. Stable identifiers and deposited data filenames are unchanged. `week` is the actual sampling week; `time_label` is its display label, and `visit_order` retains the original collection order. Repeated specimens are from the same 10 mice, not independent animals.

The expected exports comprise 851 `present` files and 24 `observed_zero_no_export` slots. Only the intended BCR or TCR chain-export slots are included; non-target and unclassified chain-report records are outside the scope of this table. No synthetic zero-row TSV files were created. `NA` denotes an unavailable or inapplicable value, not a measured zero. `assemble_chain_count_source` distinguishes explicit, inferred-zero and missing counts. Preserve strings when loading to retain those distinctions.

Native file row counts and assembly-report counts have different meanings and are stored separately. SeqKit weighted GC and FastQC GC likewise retain their original measurements. Native exports may include nonproductive sequences and repeated TRA/TRD assignments; use the documented analysis filters. MD5 values are stored once in the file records. Native TSVs themselves remain in their separate tissue archives.

Run the following command from the directory containing `reproduce.py` to regenerate the figures and tables.

```sh
python reproduce.py --omix /path/to/SLE_metadata_QC_PRJCA075761.zip --fastqc /path/to/FastQC_archive.zip --output results
```

The metadata ZIP may instead be supplied as an extracted directory. The separate FastQC archive remains necessary for per-position quality and trimmed-length data in Fig. 2; the 500 raw-file summaries here do not replace it. A completed run's prepared `results/data` can be reused with `--data-dir`.

Schema version is `2.0`. The accompanying code also supports the previous archive structure. Shared constants are verified values from this dataset; reconstruction equalities are not assumed to hold for other experiments. A JSON object containing `$ref` refers to the named dot-separated path within `study.json`. Historical provenance hashes refer to the named historical sources, not to these reorganized files. Missing BioSample/CRX/CRR accessions have not been invented or replaced by submission row IDs.
