#!/usr/bin/env python3
import argparse
from collections import defaultdict

def load_bed(path):
    intervals = set()
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                chrom, start, end, *_ = line.strip().split("\t")
                intervals.add((chrom, int(start), int(end)))
    return intervals

def match(preds, truths):
    tp = preds & truths
    fp = preds - truths
    fn = truths - preds
    return len(tp), len(fp), len(fn)

def score(tp, fp, fn):
    prec = tp / (tp + fp) if (tp + fp) else 0
    rec = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
    return prec, rec, f1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-bed", required=True)
    ap.add_argument("--pred-bed", required=True)
    args = ap.parse_args()

    truths = load_bed(args.truth_bed)
    preds = load_bed(args.pred_bed)
    tp, fp, fn = match(preds, truths)
    prec, rec, f1 = score(tp, fp, fn)

    print(f"TP\t{tp}")
    print(f"FP\t{fp}")
    print(f"FN\t{fn}")
    print(f"Precision\t{prec:.4f}")
    print(f"Recall\t{rec:.4f}")
    print(f"F1\t{f1:.4f}")

if __name__ == "__main__":
    main()
