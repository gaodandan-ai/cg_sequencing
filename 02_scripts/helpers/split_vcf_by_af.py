#!/usr/bin/env python3
from __future__ import annotations
import argparse
import gzip
from pathlib import Path

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def parse_info(info_str):
    d = {}
    if info_str == "." or not info_str:
        return d
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "1"
    return d

def safe_int(x):
    if x is None:
        return None
    try:
        return int(float(x))
    except Exception:
        return None

def pick_ao_dp(info, fmt_keys, sample_fields, alt_index=0):
    fmt = dict(zip(fmt_keys, sample_fields)) if fmt_keys and sample_fields else {}
    ao = None
    dp = None

    # FORMAT/DP
    if "DP" in fmt:
        dp = safe_int(fmt.get("DP"))

    # FORMAT/AO
    if "AO" in fmt:
        vals = fmt.get("AO", "").split(",")
        if alt_index < len(vals):
            ao = safe_int(vals[alt_index])

    # FORMAT/AD = ref,alt1,alt2...
    if (ao is None or dp is None) and "AD" in fmt:
        vals = [safe_int(x) for x in fmt.get("AD", "").split(",")]
        if len(vals) > 1 + alt_index:
            ao = vals[1 + alt_index]
        if dp is None:
            valid = [x for x in vals if x is not None]
            if valid:
                dp = sum(valid)

    # INFO fallback
    if dp is None and "DP" in info:
        dp = safe_int(info.get("DP"))
    if ao is None and "AO" in info:
        vals = info.get("AO", "").split(",")
        if alt_index < len(vals):
            ao = safe_int(vals[alt_index])

    return ao, dp

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--poly_sites", required=True)
    ap.add_argument("--high_sites", required=True)
    ap.add_argument("--poly_tsv", required=True)
    ap.add_argument("--high_tsv", required=True)
    ap.add_argument("--min_dp", type=int, default=20)
    ap.add_argument("--min_af_poly", type=float, default=0.01)
    ap.add_argument("--max_af_poly", type=float, default=0.99)
    ap.add_argument("--min_af_high", type=float, default=0.8)
    args = ap.parse_args()

    poly_sites = []
    high_sites = []
    poly_rows = []
    high_rows = []

    with open_text(Path(args.vcf)) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _id, ref, alts, qual, flt, info_str = cols[:8]
            fmt_keys = []
            sample_fields = []

            if len(cols) >= 10:
                fmt_keys = cols[8].split(":")
                sample_fields = cols[9].split(":")

            info = parse_info(info_str)

            for alt_index, alt in enumerate(alts.split(",")):
                ao, dp = pick_ao_dp(info, fmt_keys, sample_fields, alt_index)

                if dp is None or dp < args.min_dp:
                    continue

                if ao is None:
                    continue

                af = ao / dp if dp > 0 else 0.0
                site = f"{chrom}\t{pos}\n"
                row = [chrom, pos, ref, alt, dp, ao, f"{float(af):.6f}"]

                if args.min_af_poly <= af < args.max_af_poly:
                    poly_sites.append(site)
                    poly_rows.append(row)

                if af >= args.min_af_high:
                    high_sites.append(site)
                    high_rows.append(row)

    def write_unique_lines(path, lines):
        seen = set()
        with open(path, "w") as out:
            for x in lines:
                if x not in seen:
                    out.write(x)
                    seen.add(x)

    write_unique_lines(args.poly_sites, poly_sites)
    write_unique_lines(args.high_sites, high_sites)

    header = "CHROM\tPOS\tREF\tALT\tDP\tAO\tAF\n"
    with open(args.poly_tsv, "w") as out:
        out.write(header)
        for r in poly_rows:
            out.write("\t".join(map(str, r)) + "\n")

    with open(args.high_tsv, "w") as out:
        out.write(header)
        for r in high_rows:
            out.write("\t".join(map(str, r)) + "\n")

if __name__ == "__main__":
    main()
