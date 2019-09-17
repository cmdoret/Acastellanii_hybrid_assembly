# Use Mummer's dnadiff utility to compute similarity between all pairs of contigs and
# for every pair of highly similar contigs, remove the smallest

import os
import sys
import multiprocessing
from Bio import SeqIO
import argparse

# Add args genome, output_fasta, similarity threshold

# Split genome into 1 fa / contigs
# Run mummer on each possible pair
# Run mummer instances in parallel using a pool of processes
# Visualise distribution of similarity
# Write valid contigs to output fasta

