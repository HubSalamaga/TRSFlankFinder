#!/usr/bin/env python3
import argparse, random, sys, re
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def all_trimers():
    bases = "ACGT"
    return [a+b+c for a in bases for b in bases for c in bases]

def read_motifs_from_file(p):
    motifs = []
    with open(p) as f:
        for line in f:
            # accept tokens separated by whitespace or '#', ignore non-ACGT chars
            toks = re.split(r'[\s#]+', line.strip())
            for t in toks:
                tt = t.strip().upper()
                if len(tt)==3 and all(ch in "ACGT" for ch in tt):
                    motifs.append(tt)
    return motifs

def write_bed(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Make synthetic genome with planted 3-mer TRSs and truth BEDs.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genome-name", default="chrSynthetic")
    ap.add_argument("--genome-length", type=int, default=500000)
    ap.add_argument("--num-trs", type=int, default=1000)
    ap.add_argument("--min-tract", type=int, default=6)
    ap.add_argument("--max-tract", type=int, default=30)
    ap.add_argument("--min-gap", type=int, default=100, help="Minimum allowed inter-TRS gap (bp)")
    ap.add_argument("--max-gap", type=int, default=500, help="Maximum allowed inter-TRS gap (bp)")
    ap.add_argument("--motifs", default="ALL",
                    help="Comma-separated 3-mer motifs, or 'ALL' for all 64 trimers, or a file path ending with .txt to read motifs from.")
    ap.add_argument("--exclude-homopolymers", action="store_true",
                    help="If set, exclude AAA/CCC/GGG/TTT from the motif set.")
    ap.add_argument("--seed", type=int, default=1337)
    args = ap.parse_args()

    random.seed(args.seed)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Build random backbone genome (~equal base freq).
    bases = ["A","C","G","T"]
    genome = [random.choice(bases) for _ in range(args.genome_length)]

    # Resolve motif set
    motifs = None
    if args.motifs.upper() == "ALL":
        motifs = all_trimers()
    elif args.motifs.lower().endswith(".txt") and Path(args.motifs).exists():
        motifs = read_motifs_from_file(args.motifs)
        if not motifs:
            print("ERROR: No valid 3-mer motifs parsed from file.", file=sys.stderr)
            sys.exit(2)
    else:
        motifs = [m.strip().upper() for m in args.motifs.split(",") if m.strip()]
        for m in motifs:
            if len(m)!=3 or any(ch not in "ACGT" for ch in m):
                print(f"ERROR: Invalid motif '{m}'. Must be 3-mer over A/C/G/T.", file=sys.stderr)
                sys.exit(2)

    if args.exclude_homopolymers:
        motifs = [m for m in motifs if m not in ("AAA","CCC","GGG","TTT")]

    if not motifs:
        print("ERROR: motif set is empty.", file=sys.stderr)
        sys.exit(2)

    pos = 0
    trs_truth = []
    inter_truth = []

    # Place repeating pattern: TRS -> gap -> TRS -> spacer
    while len(trs_truth) < args.num_trs and pos < args.genome_length - max(args.max_gap, args.max_tract):
        motif = random.choice(motifs)
        tract_len = random.randint(args.min_tract, args.max_tract)
        tract_len -= (tract_len % 3)
        if tract_len < 6: tract_len = 6

        end_trs = pos + tract_len
        if end_trs >= args.genome_length: break
        for i in range(pos, end_trs):
            genome[i] = motif[(i - pos) % 3]
        trs_truth.append( (args.genome_name, pos, end_trs, motif, tract_len) )

        gap = random.randint(args.min_gap, args.max_gap)
        inter_start = end_trs
        inter_end   = end_trs + gap
        if inter_end >= args.genome_length: break

        # Right-bounding TRS
        motif2 = random.choice(motifs)
        tract_len2 = random.randint(args.min_tract, args.max_tract)
        tract_len2 -= (tract_len2 % 3)
        if tract_len2 < 6: tract_len2 = 6
        trs2_start = inter_end
        trs2_end   = trs2_start + tract_len2
        if trs2_end >= args.genome_length: break
        for i in range(trs2_start, trs2_end):
            genome[i] = motif2[(i - trs2_start) % 3]
        trs_truth.append( (args.genome_name, trs2_start, trs2_end, motif2, tract_len2) )

        # Inter-TRS truth
        inter_truth.append( (args.genome_name, end_trs, trs2_start, f"{motif}-{motif2}") )

        pos = trs2_end + random.randint(0, args.min_gap)  # small spacer

    # Emit outputs
    fasta_path = outdir / "synthetic_genome.fasta"
    trs_bed = outdir / "trs_truth.bed"
    inter_bed = outdir / "inter_trs_truth.bed"
    motifs_txt = outdir / "motifs_used.txt"

    record = SeqRecord(Seq("".join(genome)), id=args.genome_name, description="synthetic genome with planted TRSs")
    with open(fasta_path, "w") as f:
        SeqIO.write([record], f, "fasta")

    with open(trs_bed, "w") as f:
        for chrom, start, end, motif, tlen in trs_truth:
            f.write(f"{chrom}\t{start}\t{end}\t{motif}\t{tlen}\n")

    with open(inter_bed, "w") as f:
        for chrom, start, end, tag in inter_truth:
            if end > start:
                f.write(f"{chrom}\t{start}\t{end}\t{tag}\n")

    with open(motifs_txt, "w") as f:
        for m in sorted(set([m for m in motifs])):
            f.write(m + "\n")

    print(f"Wrote: {fasta_path}")
    print(f"Wrote: {trs_bed}")
    print(f"Wrote: {inter_bed}")
    print(f"Wrote: {motifs_txt}")
    print(f"Motif count used: {len(set(motifs))}")

if __name__ == "__main__":
    main()
