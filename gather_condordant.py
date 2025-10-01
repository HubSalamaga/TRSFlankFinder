#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def parse_metrics_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        for line in file:
            if ':' in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                data[key] = value

    return {
        "Sample": None,  # filled below
        "Total TRS in Truth": int(data.get("Total TRS in Truth", 0)),
        "Concordant Predictions": int(data.get("Concordant Predictions", 0)),
        "Concordance Rate": float(data.get("Concordance Rate", 0))
    }

def main():
    parser = argparse.ArgumentParser(description="Gather concordant TRS metrics from multiple folders into a summary CSV.")
    parser.add_argument("input_dir", help="Top-level directory to search for 'metrics_concordant.txt' files.")
    parser.add_argument("output_dir", help="Directory to write the summary CSV.")
    args = parser.parse_args()

    results = []
    for root, _, files in os.walk(args.input_dir):
        for file in files:
            if file == "metrics_concordant.txt":
                full_path = os.path.join(root, file)
                print(f"📄 Found: {full_path}")
                parsed = parse_metrics_file(full_path)

                parent_folder = os.path.basename(os.path.dirname(root))
                if parent_folder.endswith("_results"):
                    sample_name = parent_folder.replace("_results", "")
                else:
                    sample_name = parent_folder

                parsed["Sample"] = sample_name
                results.append(parsed)

    if not results:
        print("❌ No concordant metrics found.")
        return

    df = pd.DataFrame(results)
    os.makedirs(args.output_dir, exist_ok=True)
    output_csv = os.path.join(args.output_dir, "summary_concordant.csv")
    df.to_csv(output_csv, index=False)
    print(f"✅ Summary written to: {output_csv}")

if __name__ == "__main__":
    main()

