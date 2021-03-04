#!/usr/bin/env python3
# Script to excise mtDNA fragments from the assembly
# It takes a TSV file with one chromosome per line, and
# mtDNA insertion coordinates separated by tabs in the
# following format
# chrA	start1,end1	start2,end2	...
# chrB	startN,endN	...
# cmdoret, 20210304
# Usage: python excise_mtdna.py genome.fa mtDNA.tsv

import sys
from os.path import basename
from Bio import SeqIO

# Whether coordinates start at 1 in the table (as opposed to 0)
one_based = int(True)

if len(sys.argv) != 3:
    print(f"usage: python {basename(__file__)} genome.fa mtDNA.tsv")
    sys.exit(1)

mtDNA = {}
with open(sys.argv[2]) as mt:
    for line in mt:
        # Split chrom and coords and remove newline
        line = line.strip().split('\t')
        chrom = line[0]
        line = line[1:]
        # Insertions sorted by coordinates.
        # We process them in reverse order to avoid shifting coords
        insertions = line[::-1]
        mtDNA[chrom] = insertions

for scf in SeqIO.parse(sys.argv[1], 'fasta'):
    seq = str(scf.seq)
    if scf.id in mtDNA:
        for ins in mtDNA[scf.id]:
            start, end = ins.split(',')
            # Decrease coord by 1 if one based
            start = int(start) - one_based
            end = int(end) - one_based
            seq = seq[:int(start)] + seq[int(end):]
    print(">" + scf.id)
    print(seq)
