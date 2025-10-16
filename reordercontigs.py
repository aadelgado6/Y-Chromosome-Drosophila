#!/usr/bin/env python3
"""

Usage:
  reordercontigs.py -f input.fasta -o order.txt -O reordered.fasta
What it does:
  - Reads the FASTA (handles multi-line sequences)
  - Reads order.txt (one header per line, with or without leading '>')
  - Writes contigs in the exact order given in order.txt (if found)
  - Appends any contigs not listed in order.txt in their original FASTA order
  - Prints a short summary of matches / misses
"""

import argparse
import textwrap
import sys

def parse_fasta(path):
    headers = []
    seqs = {}
    cur_h = None
    cur_seq_parts = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if cur_h is not None:
                    seqs[cur_h] = ''.join(cur_seq_parts)
                cur_h = line[1:].strip()            # store header without leading '>'
                headers.append(cur_h)
                cur_seq_parts = []
            else:
                cur_seq_parts.append(line.strip())
        if cur_h is not None:
            seqs[cur_h] = ''.join(cur_seq_parts)
    return headers, seqs

def write_seq(outfh, header, seq, wrap):
    outfh.write('>' + header + '\n')
    if wrap and wrap > 0:
        for i in range(0, len(seq), wrap):
            outfh.write(seq[i:i+wrap] + '\n')
    else:
        outfh.write(seq + '\n')

def read_order_list(path):
    order = []
    with open(path, 'r') as fh:
        for ln in fh:
            s = ln.strip()
            if not s or s.startswith('#'):
                continue
            if s.startswith('>'):
                s = s[1:].strip()
            order.append(s)
    return order

def main():
    p = argparse.ArgumentParser(description="Reorder FASTA by a list and append remaining contigs.")
    p.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    p.add_argument('-o', '--order', required=True, help='Order file (one header per line, > optional)')
    p.add_argument('-O', '--output', required=True, help='Output reordered FASTA')
    p.add_argument('-w', '--wrap', type=int, default=60, help='Wrap sequence lines to N columns (default 60). Use 0 for no wrap.')
    args = p.parse_args()

    headers, seqs = parse_fasta(args.fasta)
    if not headers:
        print("No sequences found in fasta:", args.fasta, file=sys.stderr)
        sys.exit(1)

    # map first token -> list of headers (to support "contig1" matching a "contig1 something" header)
    name_to_headers = {}
    for h in headers:
        first = h.split()[0]
        name_to_headers.setdefault(first, []).append(h)

    order_list = read_order_list(args.order)

    ordered_headers = []
    used = set()
    not_found = []

    for entry in order_list:
        # try exact header match first
        if entry in seqs and entry not in used:
            ordered_headers.append(entry)
            used.add(entry)
            continue
        # try match by first token
        if entry in name_to_headers:
            # pick first header with that name that hasn't already been used
            choice = None
            for ch in name_to_headers[entry]:
                if ch not in used:
                    choice = ch
                    break
            if choice:
                ordered_headers.append(choice)
                used.add(choice)
                continue
        # not matched
        not_found.append(entry)

    # prepare output
    with open(args.output, 'w') as outfh:
        # write ordered entries first
        for h in ordered_headers:
            write_seq(outfh, h, seqs[h], args.wrap)
        # then append remaining headers in original FASTA order
        appended = 0
        for h in headers:
            if h not in used:
                write_seq(outfh, h, seqs[h], args.wrap)
                appended += 1

    # summary
    print(f"Input FASTA: {args.fasta}")
    print(f"Order list:  {args.order}")
    print(f"Output FASTA: {args.output}")
    print(f"Ordered (from order.txt) written: {len(ordered_headers)}")
    print(f"Remaining contigs appended: {appended}")
    if not_found:
        print(f"Order entries not found in FASTA ({len(not_found)}):")
        # print up to 50, otherwise truncate
        for nf in not_found[:50]:
            print("  " + nf)
        if len(not_found) > 50:
            print(f"  ... and {len(not_found)-50} more")
    else:
        print("All order.txt entries were found.")

if __name__ == '__main__':
    main()

