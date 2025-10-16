#!/usr/bin/env python3
import sys, re

if len(sys.argv) != 2:
    print("Usage: find_Ns.py assembly.fasta", file=sys.stderr)
    sys.exit(1)

fasta = sys.argv[1]
out_tsv = fasta + ".N_runs.tsv"
out_bed = fasta + ".N_runs.bed"

def parse_fasta(fp):
    name = None
    seq_lines = []
    for line in fp:
        line = line.rstrip()
        if not line: continue
        if line.startswith(">"):
            if name:
                yield name, "".join(seq_lines)
            name = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line.upper())
    if name:
        yield name, "".join(seq_lines)

with open(fasta) as fh, open(out_tsv, "w") as tsv, open(out_bed, "w") as bed:
    tsv.write("contig\tcontig_len\ttotal_Ns\tfrac_N\tn_N_runs\tn_run_positions\n")
    bed.write("#chrom\tstart\tend\tname\tN_len\n")
    for name, seq in parse_fasta(fh):
        seqlen = len(seq)
        total_Ns = seq.count("N")
        if total_Ns == 0:
            continue
        frac = total_Ns / seqlen
        runs = []
        # find runs of consecutive Ns (1-based positions)
        for m in re.finditer(r"N+", seq):
            s = m.start() + 1
            e = m.end()     # inclusive end in 1-based coords
            runs.append(f"{s}-{e}")
            # BED is 0-based half-open
            bed.write(f"{name}\t{m.start()}\t{m.end()}\t{name}_Nrun_{s}_{e}\t{m.end()-m.start()}\n")
        tsv.write(f"{name}\t{seqlen}\t{total_Ns}\t{frac:.6f}\t{len(runs)}\t{';'.join(runs)}\n")

print(f"Wrote {out_tsv} and {out_bed}")

