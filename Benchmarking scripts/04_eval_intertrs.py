#!/usr/bin/env python3
import argparse

def load_bed(path):
    regions = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                chrom, start, end, *_ = line.strip().split("\t")
                regions.append((chrom, int(start), int(end)))
    return regions

def jaccard(a, b):
    inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    union = max(a[2], b[2]) - min(a[1], b[1])
    return inter / union if union else 0

def match(preds, truths, threshold=0.5):
    matched = []
    tp = 0
    for p in preds:
        for t in truths:
            if p[0] == t[0] and jaccard(p, t) >= threshold:
                tp += 1
                matched.append((p, t))
                break
    fp = len(preds) - tp
    fn = len(truths) - tp
    return tp, fp, fn, matched

def score(tp, fp, fn):
    prec = tp / (tp + fp) if (tp + fp) else 0
    rec = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
    return prec, rec, f1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-inter-bed", required=True)
    ap.add_argument("--pred-inter-bed", required=True)
    args = ap.parse_args()

    truths = load_bed(args.truth_inter_bed)
    preds = load_bed(args.pred_inter_bed)
    tp, fp, fn, matched = match(preds, truths)
    prec, rec, f1 = score(tp, fp, fn)
    mean_jaccard = sum(jaccard(p, t) for p, t in matched) / len(matched) if matched else 0

    print(f"TP\t{tp}")
    print(f"FP\t{fp}")
    print(f"FN\t{fn}")
    print(f"Precision\t{prec:.4f}")
    print(f"Recall\t{rec:.4f}")
    print(f"F1\t{f1:.4f}")
    print(f"MeanJaccard\t{mean_jaccard:.4f}")
    print(f"Matches\t{len(matched)}")

if __name__ == "__main__":
    main()
