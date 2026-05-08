#!/usr/bin/env python3
from __future__ import annotations
import argparse
from pathlib import Path
import re
import sys

PATTERNS = [
    (re.compile(r"^(?P<sample>.+)_1\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_2\.(fastq|fq)\.gz$"), "R2"),
    (re.compile(r"^(?P<sample>.+)_R1\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_R2\.(fastq|fq)\.gz$"), "R2"),
    (re.compile(r"^(?P<sample>.+)_R1_001\.(fastq|fq)\.gz$"), "R1"),
    (re.compile(r"^(?P<sample>.+)_R2_001\.(fastq|fq)\.gz$"), "R2"),
]

def detect_pairs(indir: Path):
    samples = {}
    unmatched = []

    files = sorted(set(list(indir.glob("*.fastq.gz")) + list(indir.glob("*.fq.gz"))))

    for p in files:
        name = p.name
        matched = False
        for rx, read_tag in PATTERNS:
            m = rx.match(name)
            if m:
                sample = m.group("sample")
                samples.setdefault(sample, {})
                if read_tag in samples[sample]:
                    print(f"ERROR: duplicated {read_tag} for sample {sample}: {name}", file=sys.stderr)
                    sys.exit(1)
                samples[sample][read_tag] = str(p.resolve())
                matched = True
                break
        if not matched:
            unmatched.append(name)

    good = []
    bad = []

    for sample, d in sorted(samples.items()):
        r1 = d.get("R1")
        r2 = d.get("R2")
        if r1 and r2:
            good.append((sample, r1, r2))
        else:
            bad.append((sample, r1, r2))

    return good, bad, unmatched

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    good, bad, unmatched = detect_pairs(Path(args.indir))

    if unmatched:
        print("WARNING: unmatched FASTQ files:", file=sys.stderr)
        for x in unmatched:
            print(f"  {x}", file=sys.stderr)

    if bad:
        print("ERROR: some samples do not have complete R1/R2 pairs:", file=sys.stderr)
        for sample, r1, r2 in bad:
            print(f"  {sample}\tR1={r1}\tR2={r2}", file=sys.stderr)
        sys.exit(1)

    if not good:
        print("ERROR: no paired FASTQ files detected.", file=sys.stderr)
        sys.exit(1)

    with open(args.out, "w") as out:
        out.write("sample\tfq1\tfq2\n")
        for sample, fq1, fq2 in good:
            out.write(f"{sample}\t{fq1}\t{fq2}\n")

if __name__ == "__main__":
    main()
