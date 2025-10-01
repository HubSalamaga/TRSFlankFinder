#!/usr/bin/env python3
import argparse, re, sys

def parse_out(infile, seq_name):
    """
    Parse RepeatMasker .out format.
    We select entries with class 'Simple_repeat' and a motif in the repeat name like '(CAG)n' of length 3.
    Columns are space-delimited but positions may shift; we use a tolerant approach:
      - skip header lines until a line starting with numbers is seen after a line starting with 'SW'
      - capture fields by splitting on whitespace
    Expected fields (RepeatMasker 4.x):
      [0] SW score
      [1] perc div
      [2] perc del
      [3] perc ins
      [4] query sequence
      [5] q begin
      [6] q end
      [7] q (left)
      [8] strand (optional as 'C' for complement)
      [9] matching repeat (name)  e.g. '(CAG)n'
      [10] class/family          e.g. 'Simple_repeat'
      ... followed by repeat positions
    """
    results = []
    header_seen = False
    with open(infile) as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("SW "):  # header marker
                header_seen = True
                continue
            if not header_seen:
                continue
            toks = line.strip().split()
            if len(toks) < 11:
                continue
            # Some outputs include a 'C' complement flag shifting columns by 1
            # Detect class/family by scanning tokens for 'Simple_repeat'
            cls_idx = None
            for i, t in enumerate(toks):
                if t == "Simple_repeat":
                    cls_idx = i
                    break
            if cls_idx is None:
                continue
            # Repeat name should be token just before class/family or one/two tokens before (handles '(CAG)n')
            repname = toks[cls_idx-1]
            # Query coords usually at fixed positions after 'query sequence' which we find as the first non-numeric token
            # However safer: Get numbers that look like start<=end before class index.
            # We'll heuristically take the last two integers before cls_idx as qBegin,qEnd.
            qnums = [int(t) for t in toks[:cls_idx] if re.fullmatch(r"-?\d+", t)]
            if len(qnums) < 2:
                continue
            qBegin, qEnd = qnums[-2], qnums[-1]
            # .out is 1-based inclusive; convert to 0-based half-open
            start = min(qBegin, qEnd) - 1
            end   = max(qBegin, qEnd)
            m = re.search(r'\(([ACGT]+)\)n', repname, re.IGNORECASE)
            if not m:
                continue
            motif = m.group(1).upper()
            if len(motif) == 3:
                results.append( (seq_name, start, end, motif) )
    return results

def parse_gff(infile, default_seq_name=None):
    """
    Parse RepeatMasker GFF3.
    We look for features of type 'repeat_region' with 'Target' or 'Name' indicating Simple_repeat and try to recover motif (e.g., motif=CAG or RepeatMasker_name=(CAG)n).
    We keep only motif length 3.
    """
    results = []
    with open(infile) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = toks
            if ftype.lower() not in ("repeat_region", "simple_repeat", "region"):
                # Many RM GFFs use repeat_region for all hits
                pass
            # Require Simple_repeat in attributes if present
            if "Simple_repeat" not in attrs and "simple_repeat" not in attrs.lower():
                # Some GFFs may not label class; still try to extract motif
                pass
            # Try motif= or Repeat name like (CAG)n
            motif = None
            m1 = re.search(r'[;:]motif=([ACGT]+)', attrs, re.IGNORECASE)
            if m1:
                motif = m1.group(1).upper()
            else:
                m2 = re.search(r'[(]([ACGT]+)[)]n', attrs, re.IGNORECASE)
                if m2:
                    motif = m2.group(1).upper()
            if motif and len(motif) == 3:
                s = int(start) - 1
                e = int(end)
                results.append( (seqid or default_seq_name or ".", s, e, motif) )
    return results

def main():
    ap = argparse.ArgumentParser(description="Convert RepeatMasker .out or .gff to BED of trinucleotide Simple_repeats.")
    ap.add_argument("--infile", required=True, help=".out or .gff")
    ap.add_argument("--seq-name", default=None, help="Override sequence name (for .out)")
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    out_rows = []
    if args.infile.endswith(".gff") or args.infile.endswith(".gff3"):
        out_rows = parse_gff(args.infile, args.seq_name)
    else:
        if not args.seq_name:
            print("WARNING: --seq-name not provided; using 'seq' for .out parsing", file=sys.stderr)
        out_rows = parse_out(args.infile, args.seq_name or "seq")

    with open(args.out_bed, "w") as out:
        for chrom, s, e, motif in out_rows:
            out.write(f"{chrom}\t{s}\t{e}\t{motif}\n")

    print(f"Wrote {len(out_rows)} entries to {args.out_bed}")

if __name__ == "__main__":
    main()
