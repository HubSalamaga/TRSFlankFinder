#!/usr/bin/env python3
import argparse

def load_bed(path):
    trs = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                chrom, start, end, *rest = line.strip().split("\t")
                trs.append((chrom, int(start), int(end)))
    trs.sort()
    return trs

def find_intervals(trs, min_gap, max_gap):
    inter = []
    for i in range(len(trs) - 1):
        c1, e1, s2 = trs[i][0], trs[i][2], trs[i+1][1]
        if trs[i][0] != trs[i+1][0]:
            continue
        gap = trs[i+1][1] - trs[i][2]
        if min_gap <= gap <= max_gap:
            inter.append((c1, e1, s2))
    return inter

def write_bed(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trs-bed", required=True)
    parser.add_argument("--min-gap", type=int, default=100)
    parser.add_argument("--max-gap", type=int, default=500)
    parser.add_argument("--out-bed", required=True)
    args = parser.parse_args()

    trs = load_bed(args.trs_bed)
    inter = find_intervals(trs, args.min_gap, args.max_gap)
    write_bed(args.out_bed, inter)

if __name__ == "__main__":
    main()
