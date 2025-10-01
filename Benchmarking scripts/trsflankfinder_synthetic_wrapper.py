#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fasta_folder_path", required=True)
    parser.add_argument("--output_folder_path", required=True)
    parser.add_argument("--length_to_extract", type=int, default=100)
    parser.add_argument("--tmin", type=int, default=100)
    parser.add_argument("--tmax", type=int, default=500)
    parser.add_argument("--threshold", type=int, default=90)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--mode", choices=["combined"], default="combined")
    parser.add_argument("--keep_intermediate", action="store_true")
    args = parser.parse_args()

    os.makedirs(args.output_folder_path, exist_ok=True)

    cmd = [
        "python3", "TRSFlankFinder.py",
        "--input_fasta_folder_path", args.input_fasta_folder_path,
        "--output_folder_path", args.output_folder_path,
        "--length_to_extract", str(args.length_to_extract),
        "--tmin", str(args.tmin),
        "--tmax", str(args.tmax),
        "--threshold", str(args.threshold),
        "--threads", str(args.threads),
        "--mode", args.mode,
        "--skip_blast"
    ]

    if args.keep_intermediate:
        cmd.append("--keep_intermediate")

    print("Running TRSFlankFinder up to pre-BLAST filtering step...")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
