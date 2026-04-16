#!/usr/bin/env python3
"""
add_exons_from_cds.py

For any mRNA in a GFF3 that has CDS children but no exon children,
synthesise exon features from the CDS coordinates.

Usage:
    python add_exons_from_cds.py input.gff output.gff
"""

import sys
from collections import defaultdict


def parse_attrs(attr_str):
    attrs = {}
    for field in attr_str.strip().split(';'):
        field = field.strip()
        if '=' in field:
            k, v = field.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


def attrs_to_str(attrs):
    return ';'.join(f"{k}={v}" for k, v in attrs.items())


def process(input_path, output_path):
    records = []
    header_lines = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('#') or not line.strip():
                header_lines.append(line)
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_str = fields
            attrs = parse_attrs(attr_str)
            records.append({
                'seqid': seqid, 'source': source, 'ftype': ftype,
                'start': int(start), 'end': int(end), 'score': score,
                'strand': strand, 'phase': phase, 'attrs': attrs,
                'raw_fields': fields
            })

    # Find which mRNAs already have exons
    mrnas_with_exons = {
        r['attrs']['Parent']
        for r in records
        if r['ftype'] == 'exon' and 'Parent' in r['attrs']
    }

    # Find CDS features whose parent mRNA has no exons
    synthetic_exons = []
    exon_counter = defaultdict(int)

    for r in records:
        if r['ftype'] == 'CDS' and 'Parent' in r['attrs']:
            parent = r['attrs']['Parent']
            if parent not in mrnas_with_exons:
                exon_counter[parent] += 1
                cds_id = r['attrs'].get('ID', f"{parent}.cds{exon_counter[parent]}")
                exon_id = cds_id.replace('.cds', '.exon').replace('cds-', 'exon-')
                exon_attrs = {'ID': exon_id, 'Parent': parent}
                synthetic_exons.append({
                    'seqid': r['seqid'], 'source': r['source'],
                    'ftype': 'exon',
                    'start': r['start'], 'end': r['end'],
                    'score': r['score'], 'strand': r['strand'],
                    'phase': '.', 'attrs': exon_attrs,
                })

    print(f"[info] Synthesised {len(synthetic_exons)} exon feature(s) "
          f"across {len(set(e['attrs']['Parent'] for e in synthetic_exons))} transcript(s)")

    all_records = records + synthetic_exons
    all_records.sort(key=lambda r: (r['seqid'], r['start']))

    with open(output_path, 'w') as fout:
        fout.write("##gff-version 3\n")
        for r in all_records:
            fout.write('\t'.join([
                r['seqid'], r['source'], r['ftype'],
                str(r['start']), str(r['end']),
                r['score'], r['strand'], r['phase'],
                attrs_to_str(r['attrs'])
            ]) + '\n')

    print(f"[done] Written to {output_path}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python add_exons_from_cds.py input.gff output.gff")
        sys.exit(1)
    process(sys.argv[1], sys.argv[2])
