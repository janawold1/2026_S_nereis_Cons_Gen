#!/usr/bin/env python3
"""
filter_augustus_vep.py

Filters an AUGUSTUS GFF3 for confident gene models and reformats
for compatibility with Ensembl VEP.

Filters applied:
  - Gene score >= min_score (default 0.5)
  - Optionally require multi-exon CDS (--multi-exon)
  - Optionally require minimum hint support (--min-hints, 0.0-1.0)

Reformatting applied:
  - 'transcript' -> 'mRNA'
  - stop_codon, start_codon, intron lines dropped
  - AUGUSTUS comment blocks stripped
  - ##gff-version 3 pragma written once
  - Features sorted by seqid, start position

Usage:
    python filter_augustus_vep.py input.augustus.gff output.vep.gff

    # Stricter: score >= 0.75, must be multi-exon
    python filter_augustus_vep.py input.gff output.gff --min-score 0.75 --multi-exon

    # Also require at least 50% hint support
    python filter_augustus_vep.py input.gff output.gff --min-score 0.5 --min-hints 0.5
"""

import re
import sys
import argparse
from collections import defaultdict


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

# Regex to extract hint support percentage from AUGUSTUS comment blocks
HINT_RE = re.compile(
    r'%\s*of transcript supported by hints \(any source\):\s*(\d+)'
)


def parse_augustus_gff(path):
    """
    Parse AUGUSTUS GFF3, extracting:
      - Feature records
      - Per-transcript hint support percentages from comment blocks

    Returns:
      records         : list of record dicts
      transcript_hints: dict transcript_id -> hint_pct (0-100)
    """
    SKIP_TYPES = {'stop_codon', 'start_codon', 'intron'}
    AUGUSTUS_EXON_TYPES = {'initial', 'internal', 'terminal', 'single'}

    records = []
    transcript_hints = {}
    current_transcript = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')

            # Extract hint support from comment blocks
            if line.startswith('#'):
                # Track current transcript from comment markers
                m = re.search(r'start gene (\S+)', line)
                if m:
                    current_transcript = None  # reset; will be set when we see transcript feature

                m = HINT_RE.search(line)
                if m and current_transcript:
                    transcript_hints[current_transcript] = int(m.group(1))
                continue

            if not line.strip():
                continue

            fields = line.split('\t')
            if len(fields) != 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr_str = fields

            # Skip unwanted feature types
            if ftype in SKIP_TYPES:
                continue

            # Normalise AUGUSTUS exon subtypes
            if ftype in AUGUSTUS_EXON_TYPES:
                ftype = 'exon'

            # Normalise transcript -> mRNA
            if ftype == 'transcript':
                ftype = 'mRNA'

            # Parse attributes
            attrs = {}
            for field in attr_str.strip().split(';'):
                field = field.strip()
                if '=' in field:
                    k, v = field.split('=', 1)
                    attrs[k.strip()] = v.strip()

            # Track current transcript ID for hint extraction
            if ftype == 'mRNA' and 'ID' in attrs:
                current_transcript = attrs['ID']

            records.append({
                'seqid':  seqid,
                'source': source,
                'ftype':  ftype,
                'start':  int(start),
                'end':    int(end),
                'score':  score,
                'strand': strand,
                'phase':  phase,
                'attrs':  attrs,
            })

    return records, transcript_hints


# ---------------------------------------------------------------------------
# Gene tree helpers
# ---------------------------------------------------------------------------

def build_trees(records):
    """
    Returns:
      genes         : dict gene_id -> record
      gene_mrnas    : dict gene_id -> [mrna_ids]
      mrna_cds      : dict mrna_id -> [cds records]
      all_children  : dict parent_id -> [child records]
    """
    genes        = {}
    gene_mrnas   = defaultdict(list)
    mrna_cds     = defaultdict(list)
    all_children = defaultdict(list)

    for r in records:
        fid    = r['attrs'].get('ID')
        parent = r['attrs'].get('Parent')

        if r['ftype'] == 'gene' and fid:
            genes[fid] = r
        if r['ftype'] == 'mRNA' and fid and parent:
            gene_mrnas[parent].append(fid)
        if r['ftype'] == 'CDS' and parent:
            mrna_cds[parent].append(r)
        if parent:
            all_children[parent].append(r)

    return genes, gene_mrnas, mrna_cds, all_children


def get_all_descendants(gene_id, all_children):
    """Recursively collect all descendant records of a gene."""
    result = []
    for child in all_children.get(gene_id, []):
        result.append(child)
        cid = child['attrs'].get('ID')
        if cid:
            result.extend(get_all_descendants(cid, all_children))
    return result


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def passes_filters(gene_id, genes, gene_mrnas, mrna_cds,
                   transcript_hints, min_score, multi_exon, min_hints):
    gene = genes[gene_id]

    # Score filter
    try:
        score = float(gene['score'])
    except (ValueError, TypeError):
        score = 0.0
    if score < min_score:
        return False, f"score {score:.3f} < {min_score}"

    # Check each mRNA
    mrna_ids = gene_mrnas.get(gene_id, [])
    if not mrna_ids:
        return False, "no mRNA children"

    for mrna_id in mrna_ids:
        cds_list = mrna_cds.get(mrna_id, [])

        # Multi-exon filter
        if multi_exon and len(cds_list) < 2:
            return False, f"single-exon CDS in {mrna_id}"

        # Hint support filter
        if min_hints > 0.0:
            hint_pct = transcript_hints.get(mrna_id, 0)
            if (hint_pct / 100.0) < min_hints:
                return False, (
                    f"hint support {hint_pct}% < {int(min_hints*100)}% in {mrna_id}"
                )

    return True, None


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def record_to_line(r):
    return '\t'.join([
        r['seqid'], r['source'], r['ftype'],
        str(r['start']), str(r['end']),
        r['score'], r['strand'], r['phase'],
        ';'.join(f"{k}={v}" for k, v in r['attrs'].items())
    ])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter AUGUSTUS GFF3 for confident models and reformat for VEP."
    )
    parser.add_argument('input',  help="AUGUSTUS GFF3 input file")
    parser.add_argument('output', help="Filtered, VEP-ready GFF3 output file")
    parser.add_argument('--min-score', type=float, default=0.5,
                        help="Minimum gene score (0-1, default: 0.5)")
    parser.add_argument('--multi-exon', action='store_true',
                        help="Only keep genes with >= 2 CDS features")
    parser.add_argument('--min-hints', type=float, default=0.0,
                        help="Minimum hint support fraction (0-1, default: 0 = disabled)")
    args = parser.parse_args()

    print(f"[1/4] Parsing {args.input} ...")
    records, transcript_hints = parse_augustus_gff(args.input)
    print(f"      {len(records)} feature lines parsed")
    print(f"      {len(transcript_hints)} transcripts with hint data")

    print("[2/4] Building gene trees ...")
    genes, gene_mrnas, mrna_cds, all_children = build_trees(records)
    print(f"      {len(genes)} genes found")

    print(f"[3/4] Filtering (score >= {args.min_score}"
          + (", multi-exon only" if args.multi_exon else "")
          + (f", hints >= {int(args.min_hints*100)}%" if args.min_hints > 0 else "")
          + ") ...")

    kept_genes    = []
    dropped_genes = []

    for gid, gene in genes.items():
        ok, reason = passes_filters(
            gid, genes, gene_mrnas, mrna_cds,
            transcript_hints, args.min_score,
            args.multi_exon, args.min_hints
        )
        if ok:
            kept_genes.append(gid)
        else:
            dropped_genes.append((gid, reason))

    print(f"      Kept    : {len(kept_genes)}")
    print(f"      Dropped : {len(dropped_genes)}")

    # Score distribution of kept genes
    scores = []
    for gid in kept_genes:
        try:
            scores.append(float(genes[gid]['score']))
        except ValueError:
            pass
    if scores:
        print(f"      Score range (kept): "
              f"{min(scores):.3f} – {max(scores):.3f}, "
              f"mean {sum(scores)/len(scores):.3f}")

    print(f"[4/4] Writing {args.output} ...")

    # Collect all records for kept genes and sort
    out_records = []
    for gid in kept_genes:
        out_records.append(genes[gid])
        out_records.extend(get_all_descendants(gid, all_children))

    out_records.sort(key=lambda r: (r['seqid'], r['start']))

    with open(args.output, 'w') as fout:
        fout.write("##gff-version 3\n")
        for r in out_records:
            fout.write(record_to_line(r) + '\n')

    print(f"\n{'═' * 55}")
    print(f"  Done")
    print(f"{'═' * 55}")
    print(f"  Input genes   : {len(genes)}")
    print(f"  Output genes  : {len(kept_genes)}")
    print(f"  Retention     : {100*len(kept_genes)/max(len(genes),1):.1f}%")
    print(f"  Output file   : {args.output}")
    print(f"{'═' * 55}\n")

    if dropped_genes:
        drop_log = args.output.replace('.gff', '') + '_dropped.tsv'
        with open(drop_log, 'w') as flog:
            flog.write("gene_id\treason\n")
            for gid, reason in dropped_genes:
                flog.write(f"{gid}\t{reason}\n")
        print(f"  Dropped gene log: {drop_log}\n")


if __name__ == '__main__':
    main()
