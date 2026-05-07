#!/usr/bin/env python3
"""
filter_variants.py
==================
Filters somatic SNVs from Mutect2 output (HN002.vcf.gz) using:
  1. Mutect2 FILTER field (keep PASS only)
  2. FACETS copy number segments (flag low-purity / high-LOH regions)
  3. Basic allele frequency and depth thresholds

Usage:
    python filter_variants.py \
        --vcf HN002.vcf.gz \
        --facets_seg HG002-T2_hisens.seg \
        --facets_qc  HG002-T2.qc.txt \
        --output     HN002_filtered.vcf

Dependencies:
    pip install cyvcf2 pandas
"""

import argparse
# import cyvcf2
import sys
import pandas as pd
import gzip

# ── Argument parsing ──────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter Mutect2 VCF using FACETS copy number output"
    )
    parser.add_argument("--vcf",         required=True, help="Mutect2 VCF (can be .gz)")
    parser.add_argument("--facets_seg",  required=True, help="FACETS hisens .seg file")
    parser.add_argument("--facets_qc",   required=True, help="FACETS QC file (.qc.txt)")
    parser.add_argument("--output",      required=True, help="Output filtered VCF path")
    parser.add_argument("--min_af",      type=float, default=0.05,
                        help="Minimum tumor allele frequency (default: 0.05)")
    parser.add_argument("--min_depth",   type=int,   default=10,
                        help="Minimum read depth in tumor (default: 10)")
    parser.add_argument("--max_tcn",     type=int,   default=8,
                        help="Max total copy number allowed (default: 8, flags amplifications)")
    return parser.parse_args()

# ── Load FACETS QC ────────────────────────────────────────────────────────────
def load_facets_qc(qc_path):
    """
    Parse the FACETS QC file to extract purity and ploidy estimates.
    Returns a dict with keys: purity, ploidy, dipLogR
    """
    qc = {}
    with open(qc_path) as f:
        lines = f.readlines()
    # QC file is tab-separated with a header row
    if len(lines) >= 2:
        headers = lines[0].strip().split("\t")
        values  = lines[1].strip().split("\t")
        qc = dict(zip(headers, values))
    purity = float(qc.get("purity", 0))
    ploidy = float(qc.get("ploidy", 2))
    print(f"[INFO] FACETS QC - Purity: {purity:.2f}, Ploidy: {ploidy:.2f}")
    if purity < 0.3:
        print("[WARN] Low tumor purity (<0.30). Many true somatic variants may be missed or unreliable.")
    return purity, ploidy

# ── Load FACETS segments ──────────────────────────────────────────────────────
def load_facets_segments(seg_path):
    """
    Load the FACETS hisens segmentation file.
    Columns include: chrom, loc.start, loc.end, tcn.em (total CN), lcn.em (minor CN), cf.em (cell fraction)
    Returns a pandas DataFrame.
    """
    seg = pd.read_csv(seg_path, sep="\t")
    # Normalize chromosome naming (ensure 'chr' prefix is absent for matching)
    seg["chrom"] = seg["chrom"].astype(str).str.replace("chr", "", regex=False)
    print(f"[INFO] Loaded {len(seg)} FACETS segments from {seg_path}")
    return seg

# ── Look up FACETS segment for a variant ─────────────────────────────────────
def get_facets_segment(chrom, pos, seg_df):
    """
    Given a chromosome and position, find the overlapping FACETS segment.
    Returns the segment row as a pandas Series, or None if not found.
    """
    chrom = str(chrom).replace("chr", "")
    match = seg_df[
        (seg_df["chrom"] == chrom) &
        (seg_df["loc.start"] <= pos) &
        (seg_df["loc.end"]   >= pos)
    ]
    if match.empty:
        return None
    return match.iloc[0]

# ── Parse a VCF line ──────────────────────────────────────────────────────────
def parse_vcf_line(line):
    """
    Parse a single VCF data line into its fields.
    Returns a dict with keys: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples
    """
    fields = line.strip().split("\t")
    if len(fields) < 9:
        return None
    return {
        "CHROM":   fields[0],
        "POS":     int(fields[1]),
        "ID":      fields[2],
        "REF":     fields[3],
        "ALT":     fields[4],
        "QUAL":    fields[5],
        "FILTER":  fields[6],
        "INFO":    fields[7],
        "FORMAT":  fields[8],
        "samples": fields[9:]   # [normal_sample, tumor_sample] by Mutect2 convention
    }

# ── Extract tumor allele frequency and depth ──────────────────────────────────
def get_tumor_af_depth(vcf_record):
    """
    Extract tumor allele frequency (AF) and depth (DP) from the tumor sample column.
    Mutect2 FORMAT fields include: GT:AD:AF:DP:F1R2:F2R1:FAD:SB
    The TUMOR sample is the SECOND sample column in Mutect2 output.
    Returns (af, depth) as floats, or (None, None) on failure.
    """
    fmt_keys = vcf_record["FORMAT"].split(":")
    if len(vcf_record["samples"]) < 2:
        return None, None
    tumor_vals = vcf_record["samples"][1].split(":")   # index 1 = tumor
    fmt = dict(zip(fmt_keys, tumor_vals))

    try:
        af    = float(fmt.get("AF", 0))
        depth = int(fmt.get("DP", 0))
    except (ValueError, TypeError):
        return None, None
    return af, depth

# ── Main filtering logic ──────────────────────────────────────────────────────
def filter_vcf(args):
    purity, ploidy = load_facets_qc(args.facets_qc)
    seg_df         = load_facets_segments(args.facets_seg)

    # Counters for reporting
    total = passed_filter = failed_mutect = failed_af = failed_depth = 0
    failed_cn = failed_loh = 0

    # Open input VCF (supports .gz)
    opener = gzip.open if args.vcf.endswith(".gz") else open
    mode   = "rt"

    with opener(args.vcf, mode) as vcf_in, open(args.output, "w") as vcf_out:
        for line in vcf_in:

            # ── Pass header lines through unchanged ──────────────────────────
            if line.startswith("#"):
                vcf_out.write(line)
                # Add a custom FILTER header entry for our FACETS-based filters
                if line.startswith("#CHROM"):
                    vcf_out.write('##FILTER=<ID=high_CN,Description="Total copy number > max_tcn threshold from FACETS">\n')
                    vcf_out.write('##FILTER=<ID=LOH_region,Description="Loss of heterozygosity region (lcn.em == 0) from FACETS">\n')
                    vcf_out.write('##FILTER=<ID=low_AF,Description="Tumor allele frequency below minimum threshold">\n')
                    vcf_out.write('##FILTER=<ID=low_depth,Description="Tumor read depth below minimum threshold">\n')
                continue

            total += 1
            record = parse_vcf_line(line)
            if record is None:
                continue

            filter_flags = []

            # ── Filter 1: Mutect2 PASS filter ───────────────────────────────
            # Only keep variants that passed all of Mutect2's internal filters
            # (strand bias, base quality, read orientation, contamination, etc.)
            if record["FILTER"] not in ("PASS", "."):
                failed_mutect += 1
                continue   # hard exclude — no point in further evaluation

            # ── Filter 2: Tumor allele frequency ────────────────────────────
            # Low AF variants are likely sequencing noise or germline contamination
            af, depth = get_tumor_af_depth(record)
            if af is not None and af < args.min_af:
                failed_af += 1
                filter_flags.append("low_AF")

            # ── Filter 3: Minimum read depth ────────────────────────────────
            # Low depth means insufficient evidence to call a somatic variant
            if depth is not None and depth < args.min_depth:
                failed_depth += 1
                filter_flags.append("low_depth")

            # ── Filter 4: FACETS copy number integration ─────────────────────
            seg = get_facets_segment(record["CHROM"], record["POS"], seg_df)
            if seg is not None:
                tcn = seg.get("tcn.em")   # total copy number
                lcn = seg.get("lcn.em")   # minor (lesser) copy number

                # Flag regions with extreme amplification — these regions have
                # complex allele structures that make AF interpretation unreliable
                if pd.notna(tcn) and tcn > args.max_tcn:
                    failed_cn += 1
                    filter_flags.append("high_CN")

                # Flag LOH regions (lcn == 0 means all copies are from one allele)
                # Somatic variants in LOH regions can appear at unexpected AFs
                if pd.notna(lcn) and lcn == 0:
                    failed_loh += 1
                    filter_flags.append("LOH_region")

            # ── Write output ─────────────────────────────────────────────────
            if filter_flags:
                # Update the FILTER field with our new flags but still write the variant
                # (soft filter — lets downstream tools decide)
                fields      = line.strip().split("\t")
                fields[6]   = ";".join(filter_flags)
                vcf_out.write("\t".join(fields) + "\n")
            else:
                # Variant passed all filters — write as PASS
                passed_filter += 1
                vcf_out.write(line)

    # ── Summary report ───────────────────────────────────────────────────────
    print("\n===== Filtering Summary =====")
    print(f"  Total variants evaluated : {total}")
    print(f"  Failed Mutect2 FILTER    : {failed_mutect}")
    print(f"  Failed min AF ({args.min_af})    : {failed_af}")
    print(f"  Failed min depth ({args.min_depth})   : {failed_depth}")
    print(f"  Failed high copy number  : {failed_cn}")
    print(f"  Failed LOH region        : {failed_loh}")
    print(f"  Final PASS variants      : {passed_filter}")
    print(f"  Output written to        : {args.output}")
    print("=============================\n")

# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    args = parse_args()
    filter_vcf(args)
