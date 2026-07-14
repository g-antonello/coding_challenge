#!/usr/bin/env python3
"""
Scan a samtools mpileup TSV across every position it contains and report
only the positions where a variant is called, based on two fixed rules:

  1. Coverage at the position must be >= --min-coverage (default 10X)
  2. The non-reference (mutant) read fraction must be >= --min-variant-freq
     (default 0.05)

For each qualifying position, the dominant alternate allele (the
non-reference base seen most often) is reported alongside its count and
frequency. Depth and allele counts are computed by parsing the mpileup
read-base string directly (not column 4), so indel and start/end-of-read
markers don't inflate counts. Deletions ('*') are excluded from both the
ref and alt counts, since they aren't informative for a SNP call.

By default every position in the pileup is scanned. Use --positions to
restrict reporting to a specific comma-separated list of 1-based
positions (e.g. --positions 2268,2269) if you only care about known
sites, without having to re-run mpileup with a narrower region.
"""
import argparse
import csv


BASES = "ACGTN"


def parse_pileup_bases(bases_field, ref_base):
    """Return (ref_count, alt_base_counts) from an mpileup read-bases field.

    alt_base_counts is a dict of {base: count} for A/C/G/T/N bases that
    differ from the reference (case-insensitive; strand is not tracked).
    """
    ref_count = 0
    alt_base_counts = {b: 0 for b in BASES}
    i = 0
    n = len(bases_field)
    while i < n:
        c = bases_field[i]
        if c == '^':
            # start-of-read marker, followed by one mapping-quality char
            i += 2
            continue
        if c == '$':
            # end-of-read marker
            i += 1
            continue
        if c in '.,':
            ref_count += 1
            i += 1
            continue
        if c == '*':
            # deletion placeholder at this position - not ref, not a SNP alt
            i += 1
            continue
        if c in '+-':
            # indel: c, then digits giving length, then that many bases to skip
            i += 1
            digits = ''
            while i < n and bases_field[i].isdigit():
                digits += bases_field[i]
                i += 1
            i += int(digits)
            continue
        cu = c.upper()
        if cu in BASES:
            alt_base_counts[cu] += 1
            i += 1
            continue
        # anything else (shouldn't normally occur): skip one char
        i += 1
    return ref_count, alt_base_counts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pileup", required=True, help="samtools mpileup output (TSV)")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--min-coverage", type=int, default=10)
    parser.add_argument("--min-variant-freq", type=float, default=0.05)
    parser.add_argument("--positions", default=None,
                         help="Optional comma-separated list of 1-based positions "
                              "to restrict reporting to (default: scan everything)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    positions_filter = None
    if args.positions:
        positions_filter = {int(p) for p in args.positions.split(",")}

    rows = []
    with open(args.pileup) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            contig, pos, ref_base = fields[0], int(fields[1]), fields[2]

            if positions_filter is not None and pos not in positions_filter:
                continue

            bases_field = fields[4] if len(fields) > 4 else ""
            ref_count, alt_base_counts = parse_pileup_bases(bases_field, ref_base)

            total_alt_count = sum(alt_base_counts.values())
            depth = ref_count + total_alt_count
            alt_freq = (total_alt_count / depth) if depth > 0 else 0.0

            coverage_pass = depth >= args.min_coverage
            mutation_called = coverage_pass and (alt_freq >= args.min_variant_freq)

            if not mutation_called:
                continue

            # dominant alternate allele = most frequently observed non-ref base
            alt_allele, alt_allele_count = max(
                alt_base_counts.items(), key=lambda kv: kv[1])
            alt_allele_freq = (alt_allele_count / depth) if depth > 0 else 0.0

            rows.append({
                "sample_id": args.sample_id,
                "contig": contig,
                "position": pos,
                "ref_base": ref_base,
                "depth": depth,
                "ref_count": ref_count,
                "alt_allele": alt_allele,
                "alt_allele_count": alt_allele_count,
                "alt_allele_freq": round(alt_allele_freq, 4),
                "total_alt_count": total_alt_count,
                "total_alt_freq": round(alt_freq, 4),
            })

    fieldnames = ["sample_id", "contig", "position", "ref_base", "depth",
                  "ref_count", "alt_allele", "alt_allele_count", "alt_allele_freq",
                  "total_alt_count", "total_alt_freq"]
    with open(args.output, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()