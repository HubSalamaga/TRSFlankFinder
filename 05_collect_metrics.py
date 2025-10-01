#!/usr/bin/env python3
import argparse
import os
import csv
from pathlib import Path

def parse_metrics_file(path):
    metrics = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) == 2:
                key, value = parts
                try:
                    metrics[key.strip()] = float(value)
                except ValueError:
                    metrics[key.strip()] = value
    return metrics

def main():
    ap = argparse.ArgumentParser(description="Collect benchmark metrics from TRSFlankFinder and RepeatMasker.")
    ap.add_argument("--rootdir", required=True, help="Directory containing genome subfolders and summary.csv")
    ap.add_argument("--out", required=True, help="Output summary CSV path")
    args = ap.parse_args()

    root = Path(args.rootdir)
    summary_path = root / "summary.csv"
    if not summary_path.exists():
        print(f"ERROR: summary.csv not found in {root}")
        return

    with open(summary_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    for row in rows:
        gname = row["genome"]
        gdir = root / gname

        files = {
            "TRS_RM": gdir / "metrics_trs.txt",
            "INTER_RM": gdir / "metrics_intertrs.txt",
            "TRS_TRSF": gdir / "metrics_trs_TRSF.txt",
            "INTER_TRSF": gdir / "metrics_intertrs_TRSF.txt"
        }

        for prefix, path in files.items():
            if path.exists():
                parsed = parse_metrics_file(path)
                for k, v in parsed.items():
                    row[f"{prefix}_{k}"] = v

    out_fields = sorted({k for row in rows for k in row})
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[Done] Collected metrics written to {args.out}")

if __name__ == "__main__":
    main()
