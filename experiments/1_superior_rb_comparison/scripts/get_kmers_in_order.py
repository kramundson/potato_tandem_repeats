#!/cato2pool/backup-share/kramundson/software/miniconda3/bin/python
# get_kmers_in_order.py

"""
USAGE: program.py -f fasta -k 31 -o <output>

Takes a FASTA formatted repat consensus sequence as input, then returns k-mers in the
order they appear in the sequence (forward strand)
"""

import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", dest="fasta")
parser.add_argument("-k", dest="k", type=int, default=1)
parser.add_argument("-o", dest="out", type=str)
args = parser.parse_args()

k = args.k
o = open(args.out, 'w')

# Jellyfish automatically reports the lexicographically shorter of the 2 mers.
# This script does the same. In cases where the reverse complement mer is kept, the 
# starting coordinate reported for each mer is still based on the original sequence.

forward_write = 0
reverse_write = 0

with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for i in range(0, len(record.seq)-k+1, 1):
            fwd = str(record.seq[i:i+k])
            rev = str(record.seq[i:i+k].reverse_complement())
            if fwd < rev:
                o.write("\t".join([str(i), fwd, "\n"]))
            else:
                o.write("\t".join([str(i), rev, "\n"]))
