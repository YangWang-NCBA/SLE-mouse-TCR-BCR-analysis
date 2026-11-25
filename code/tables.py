"""Generate the publication tables and the paired read counts for Fig. 2b."""

import json
from pathlib import Path
import sys

import pandas as pd


CHAINS = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']
TISSUES = [
    'blood', 'spleen', 'lymph_node', 'thymus', 'brain', 'kidney',
    'liver', 'lung', 'heart', 'intestine', 'skin',
]
GROUPS = {'untreated_control': 'Control', 'IMQ_treated': 'IMQ'}


def write_tsv(data, path):
    data.to_csv(path, sep='\t', index=False, lineterminator='\n')


def paired_read_counts(data, output):
    path = data / 'fastqc/FastQC_file_summary.tsv'
    if not path.is_file():
        path = path.with_suffix('.tsv.gz')
    summary = pd.read_csv(path, sep='\t', float_precision='round_trip')
    raw = summary.loc[summary.stage.eq('raw')]
    pairs = raw.pivot(index='analysis_library_id', columns='read_direction', values='read_count')
    pairs = pairs.rename(columns={'R1': 'R1_read_count', 'R2': 'R2_read_count'}).reset_index()
    fields = ['SampleID', 'analysis_library_id', 'mouse_id', 'group_public', 'tissue', 'receptor']
    pairs = raw.loc[raw.read_direction.eq('R1'), fields].merge(
        pairs, on='analysis_library_id', validate='one_to_one')
    trimmed = summary.loc[summary.stage.eq('trimmed')].pivot(
        index='analysis_library_id', columns='read_direction', values='read_count')
    if not (pairs.R1_read_count.eq(pairs.R2_read_count).all() and trimmed.R1.eq(trimmed.R2).all()):
        raise ValueError('Paired FASTQ files have different read counts.')
    pairs = pairs.merge(trimmed[['R1']].rename(columns={'R1': 'trimmed_read_pairs'}),
                        left_on='analysis_library_id', right_index=True, validate='one_to_one')
    pairs['retained_pair_fraction'] = pairs.trimmed_read_pairs / pairs.R1_read_count
    folder = output / 'fig/Fig2_sequencing_quality/source_data'
    folder.mkdir(parents=True, exist_ok=True)
    write_tsv(pairs, folder / 'Fig2b_R1_R2_read_count_agreement.tsv')
    return pairs


def main(root, release):
    paths = {
        'library_qc': root / 'analysis/mixcr_library_qc.tsv',
        'native_counts': root / 'analysis/native_counts.tsv',
        'chain_status': root / 'metadata/chains.tsv',
    }
    qc, native, status = [pd.read_csv(p, sep='\t') for p in paths.values()]
    reads = paired_read_counts(root, release)
    libs = qc.merge(
        reads[['analysis_library_id', 'SampleID', 'R1_read_count', 'R2_read_count']],
        on=['analysis_library_id', 'SampleID'], how='left', validate='one_to_one',
    ).merge(native, on=['analysis_library_id', 'SampleID'], how='left', validate='one_to_one')
    identity = [
        'SampleID', 'analysis_library_id', 'mouse_id',
        'group_public', 'tissue', 'receptor',
    ]
    chain_detail = libs.melt(
        id_vars=identity, value_vars=CHAINS,
        var_name='chain', value_name='assembled_clonotypes',
    )
    is_bcr = chain_detail.receptor.eq('BCR') & chain_detail.chain.str.startswith('IG')
    is_tcr = chain_detail.receptor.eq('TCR') & chain_detail.chain.str.startswith('TR')
    chain_detail = chain_detail.loc[is_bcr | is_tcr].merge(
        status[['analysis_library_id', 'chain', 'chain_file_status']],
        on=['analysis_library_id', 'chain'], how='left', validate='one_to_one',
    )

    if libs[['R1_read_count', 'umi_output_diversity', *CHAINS]].isna().any().any():
        raise ValueError('A library lacks read, UMI or clonotype counts.')
    if chain_detail.chain_file_status.isna().any():
        raise ValueError('Missing specimen–chain export status.')
    absent = status.loc[status.chain_file_status.ne('present')]
    recorded_zero = (absent.chain_file_status.eq('observed_zero_no_export')
                     & absent.assemble_chain_total.eq(0)
                     & absent.assemble_chain_count_closure_verified.eq('yes'))
    absent_counts = chain_detail.loc[chain_detail.chain_file_status.ne('present'), 'assembled_clonotypes']
    if not (recorded_zero.all() and absent_counts.eq(0).all()):
        raise ValueError('An absent export lacks a recorded zero count.')

    out = release / 'table/publication'
    source = release / 'table/source_data/Table4_tissue_totals'
    out.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    detail = libs[identity + ['R1_read_count', 'R2_read_count', 'umi_output_diversity']]
    detail = detail.rename(columns={
        'R1_read_count': 'raw_read_pairs', 'umi_output_diversity': 'refined_umis',
    })
    write_tsv(detail, source / 'library_totals.tsv')
    write_tsv(chain_detail, source / 'assembled_clonotypes_by_specimen_chain.tsv')
    labels = {tissue: tissue.capitalize() for tissue in TISSUES}
    labels.update(blood='Peripheral blood', lymph_node='Lymph nodes')
    rows = []
    for tissue in TISSUES:
        for group, label in GROUPS.items():
            selected = libs.loc[libs.tissue.eq(tissue) & libs.group_public.eq(group)]
            chain_totals = chain_detail.loc[
                chain_detail.tissue.eq(tissue) & chain_detail.group_public.eq(group)
            ].groupby('chain').assembled_clonotypes.sum()
            rows.append({
                'Specimen source': labels[tissue], 'Group': label,
                'Specimens': selected.SampleID.nunique(),
                'Raw read pairs': int(selected.R1_read_count.sum()),
                'Refined UMIs': int(selected.umi_output_diversity.sum()),
                **{chain: int(chain_totals[chain]) for chain in CHAINS},
            })
    result = pd.DataFrame(rows)
    numeric = list(result.columns[2:])
    for label in GROUPS.values():
        subtotal = result.loc[result.Group.eq(label)]
        result.loc[len(result)] = {
            'Specimen source': 'Total', 'Group': label,
            **{column: int(subtotal[column].sum()) for column in numeric},
        }
    write_tsv(result, out / 'Table4_tissue_read_umi_clonotype_totals.tsv')

    tables = json.loads((root / 'tables_1_to_3.json').read_text(encoding='utf-8'))
    for item in tables:
        write_tsv(
            pd.DataFrame(item['rows'][1:], columns=item['rows'][0]),
            out / f'Table{item["number"]}.tsv',
        )
    display = result.copy()
    for column in numeric:
        display[column] = display[column].map(lambda value: f'{int(value):,}')
    note = (
        'Control denotes untreated controls and IMQ denotes imiquimod treatment, '
        'with five mice per group. Samples are biological specimens; repeated '
        'blood collections are counted separately. Control blood was collected '
        'at weeks 0 and 8, and IMQ blood at weeks 0, 4 and 8. Read pairs and refined '
        'UMIs sum the TCR and BCR libraries. Refined UMIs are output UMI diversity '
        'from refineTagsAndSort. Clonotypes are summarised under IGH, IGK, IGL, TRA, '
        'TRB, TRD and TRG; the five heavy-chain isotypes are combined under IGH. '
        'Native cloneId values with a single-chain assignment are counted once '
        'per library across exports, including productive and nonproductive '
        'records. Counts are summed across libraries without further V–J–CDR3 '
        'merging. A further 4,400 records have an empty chains field and one has '
        'a TRAD assignment without a unique TRA or TRD assignment. These are '
        'retained in the source data and excluded from the table counts.'
    )
    tables.append({
        'number': 4,
        'title': 'Table 4 Sequencing yield and MiXCR clonotypes by specimen source and group',
        'rows': [list(display.columns)] + display.values.tolist(), 'note': note,
    })
    (out / 'publication_tables.json').write_text(
        json.dumps(tables, ensure_ascii=False, indent=2), encoding='utf-8',
    )
    lines = []
    for item in tables:
        rows = item['rows']
        lines += ['## ' + item['title'], '', '| ' + ' | '.join(rows[0]) + ' |',
                  '| ' + ' | '.join(['---'] * len(rows[0])) + ' |']
        lines += [
            '| ' + ' | '.join(str(x).replace('\n', '<br>') for x in row) + ' |'
            for row in rows[1:]
        ]
        lines += ['', item['note'], '']
    (out / 'Tables1-4.md').write_text('\n'.join(lines), encoding='utf-8')
    print(result.to_string(index=False))


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
