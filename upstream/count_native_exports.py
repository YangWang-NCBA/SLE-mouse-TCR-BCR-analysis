#!/usr/bin/env python3
"""Count each native cloneId once per library and list ambiguous chain assignments."""
from pathlib import Path
from collections import Counter, defaultdict
import argparse, csv, gzip
from sequence_metrics import sample_id_from_library

CHAINS=('IGH','IGK','IGL','TRA','TRB','TRD','TRG')

def write_tsv(path,rows,fields=None):
    fields=fields or list(rows[0])
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', encoding='utf-8', newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)

def count_library(library,paths):
    sample_id = sample_id_from_library(library)
    seen = {}
    counts = Counter()
    duplicates = 0
    unknown = []
    for path in sorted(paths):
        file_ids = set()
        with path.open('rb',buffering=8*1024*1024) as handle:
            header_bytes=handle.readline()
            header=header_bytes.decode('utf-8-sig').rstrip('\r\n').split('\t')
            idx={key:header.index(key) for key in ('cloneId','chains')}
            for line in handle:
                values=line.rstrip(b'\r\n').split(b'\t')
                if values==[b'']:continue
                if len(values)!=len(header):raise ValueError(f'Malformed row in {path.name}')
                cid=int(values[idx['cloneId']]);chain=values[idx['chains']].decode('utf-8')
                if cid in file_ids:raise ValueError(f'Duplicate cloneId within {path.name}: {cid}')
                file_ids.add(cid)
                if cid in seen:
                    duplicates+=1
                    if seen[cid]!=chain:raise ValueError(f'Conflicting chain assignments in {library}: {cid}')
                else:
                    seen[cid]=chain;counts[chain]+=1
                    if chain not in CHAINS:unknown.append({'SampleID':sample_id,'analysis_library_id':library,'cloneId':cid,'chains':chain})
    result = {
        'SampleID': sample_id,
        'analysis_library_id': library,
        'export_files': len(paths),
        'distinct_assembled_clone_ids': len(seen),
        'duplicate_rows_across_exports': duplicates,
        **{chain: counts[chain] for chain in CHAINS},
        'unassigned_chain': sum(n for chain, n in counts.items() if chain not in CHAINS),
    }
    return result,unknown

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tsv-root',required=True,type=Path)
    parser.add_argument('--release-dir',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    if args.output_dir.exists():raise ValueError('Use a new output directory')
    metrics=args.release_dir/'data/analysis/repertoire_chain_metrics.tsv'
    handle=metrics.open(encoding='utf-8-sig',newline='') if metrics.is_file() else gzip.open(metrics.with_suffix('.tsv.gz'),'rt',encoding='utf-8-sig',newline='')
    with handle:
        expected_rows=list(csv.DictReader(handle,delimiter='\t'))
    expected={r['tsv_filename']:r for r in expected_rows}
    if len(expected)!=len(expected_rows):raise ValueError('Duplicate expected filename')
    paths={}
    for p in args.tsv_root.rglob('*.tsv'):
        if p.name in expected:
            if p.name in paths:raise ValueError(f'Duplicate input basename: {p.name}')
            paths[p.name]=p
    if set(paths)!=set(expected):raise ValueError('Missing expected native TSV inputs')
    libraries=defaultdict(list)
    for name,p in paths.items():libraries[expected[name]['analysis_library_id']].append(p)
    out=args.output_dir.resolve()
    if out==args.tsv_root.resolve() or args.tsv_root.resolve() in out.parents:raise ValueError('Output must be outside the native-input tree')
    out.mkdir(parents=True)
    library_rows = []
    unknown = []
    for library,files in sorted(libraries.items()):
        row, other = count_library(library, files)
        library_rows.append(row)
        unknown.extend(other)
        print(f'{len(library_rows)}/{len(libraries)} libraries complete',flush=True)
    write_tsv(out/'native_assembled_counts_by_library.tsv',library_rows)
    write_tsv(out/'unassigned_native_clonotypes.tsv.gz',unknown,fields=['SampleID','analysis_library_id','cloneId','chains'])
if __name__=='__main__':main()
