"""Extract the three Figure 2 input tables from the deposited FastQC reports."""

import argparse
from contextlib import ExitStack
import csv
import io
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import pandas as pd

from sequence_metrics import parse_fastqc_modules, parse_table, sample_id_from_library


LABELS = ['SampleID', 'analysis_library_id', 'mouse_id', 'group_public', 'tissue', 'receptor']
FILE = LABELS + ['stage', 'read_direction', 'filename']
QUALITY = ['mean_quality', 'median_quality', 'lower_quartile', 'upper_quartile',
           'percentile_10', 'percentile_90']


def report_index(libraries):
    """Match raw and Trim Galore mate filenames to the deposited library IDs."""
    for library in libraries:
        fields = {key: library[key] for key in LABELS if key != 'SampleID'}
        fields['SampleID'] = sample_id_from_library(library['analysis_library_id'])
        for stage in ('raw', 'trimmed'):
            for direction in ('R1', 'R2'):
                name = library[direction + '_filename']
                if stage == 'trimmed':
                    name = re.sub(r'(\.f(?:ast)?q\.gz)$', r'_val_' + direction[1] + r'\1', name)
                yield {**fields, 'stage': stage, 'read_direction': direction, 'filename': name}


def extract_tables(source, libraries, output):
    with ExitStack() as files:
        def table(name, columns):
            handle = files.enter_context((output / name).open('w', encoding='utf-8', newline=''))
            writer = csv.DictWriter(handle, columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            return writer

        summary = table('FastQC_file_summary.tsv', FILE + ['read_count', 'sequence_length', 'gc_percent'])
        per_base = table('raw_per_base_metrics.tsv', FILE + ['base_bin', 'base_start', 'base_end',
                         'base_midpoint', *QUALITY, 'n_percent'])
        lengths = table('Fig2d_trimmed_read_length_distribution.tsv',
                        ['SampleID', 'length_midpoint', 'read_count', 'total_reads', 'percent_per_bp',
                         *LABELS[1:], 'read_direction', 'filename'])
        archive = files.enter_context(ZipFile(source)) if source.is_file() else None
        base = source / 'fastqc' if (source / 'fastqc').is_dir() else source
        for record in report_index(libraries):
            sample = {key: record[key] for key in FILE}
            stage, library, filename = (sample[key] for key in ('stage', 'analysis_library_id', 'filename'))
            if stage not in ('raw', 'trimmed') or not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', x)
                                                          for x in (library, filename)):
                raise ValueError(f'Invalid FastQC index entry: {library}, {filename}, {stage}')
            stem = re.sub(r'\.f(?:ast)?q\.gz$', '', filename)
            path = Path(stage, library, stem + '_fastqc.zip')
            raw = archive.read('fastqc/' + path.as_posix()) if archive else (base / path).read_bytes()
            with ZipFile(io.BytesIO(raw)) as report:
                member, = (name for name in report.namelist() if name.endswith('/fastqc_data.txt'))
                text = report.read(member).decode('utf-8-sig')
            modules = {name: parse_table(lines)[1] for name, _, lines in parse_fastqc_modules(text)}
            basic = dict(modules['Basic Statistics'])
            if basic['Filename'] != filename:
                raise ValueError(f'{path} describes {basic["Filename"]}, expected {filename}')
            total = int(basic['Total Sequences'])
            summary.writerow({**sample, 'read_count': total, 'sequence_length': basic['Sequence length'],
                              'gc_percent': float(basic['%GC'])})
            if stage == 'raw':
                n_content = dict(modules['Per base N content'])
                for position, *phred in modules['Per base sequence quality']:
                    pos = int(position)
                    per_base.writerow({**sample, 'base_bin': pos, 'base_start': pos, 'base_end': pos,
                                       'base_midpoint': float(pos), **dict(zip(QUALITY, map(float, phred), strict=True)),
                                       'n_percent': float(n_content[position])})
            else:
                counts = {int(bp): int(float(n)) for bp, n in modules['Sequence Length Distribution']}
                if total <= 0 or sum(counts.get(bp, 0) for bp in range(35, 301)) != total:
                    raise ValueError(f'{filename} has no reads or includes lengths outside 35–300 bp')
                labels = {key: sample[key] for key in LABELS + ['read_direction', 'filename']}
                for bp in range(35, 301):
                    n = counts.get(bp, 0)
                    lengths.writerow({**labels, 'length_midpoint': bp, 'read_count': n,
                                      'total_reads': total, 'percent_per_bp': 100 * n / total})

    # Preserve the publication tables' decimal representation and library ordering.
    path = output / 'raw_per_base_metrics.tsv'
    per_base = pd.read_csv(path, sep='\t')
    per_base.to_csv(path, sep='\t', index=False, lineterminator='\n')
    path = output / 'Fig2d_trimmed_read_length_distribution.tsv'
    lengths = pd.read_csv(path, sep='\t', dtype={'percent_per_bp': str})
    lengths = lengths.sort_values(['analysis_library_id', 'read_direction', 'length_midpoint'])
    lengths.to_csv(path, sep='\t', index=False, lineterminator='\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastqc-source', required=True, type=Path, help='OMIX report ZIP or extracted directory')
    parser.add_argument('--output', required=True, type=Path, help='New directory for the three TSV files')
    parser.add_argument('--libraries', required=True, type=Path, help='Prepared data/metadata/libraries.tsv')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Choose a new output directory.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=args.output.parent) as temp:
        tables = Path(temp) / 'tables'
        tables.mkdir()
        with args.libraries.open(encoding='utf-8-sig', newline='') as handle:
            extract_tables(args.fastqc_source, csv.DictReader(handle, delimiter='\t'), tables)
        tables.rename(args.output)


if __name__ == '__main__':
    main()
