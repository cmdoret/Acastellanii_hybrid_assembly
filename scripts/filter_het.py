# This script uses minimap2 python binding to filter redundant haplotypes out from a genome
# assembly. It does so by finding pairs of contigs with high identitity and overlap
# Inspired by the reduction step of redundans: https://github.com/Gabaldonlab/redundans
# cmdoret, 20190924

import mappy as mp
import itertools as it
import click
from Bio import SeqIO
import sys




@click.command()
@click.option(
    "-I",
    "--min-identity",
    help="Minimum (blast-like) identity required to consider contigs homozygous",
    default=0.51,
    show_default=True,
)
@click.option(
    "-O",
    "--min-overlap",
    help="Minimum overlap required to consider contigs homozygous. Defined as: length(alignment) / length(c) where c is the smallest contig of the aligned pair.",
    default=0.8,
    show_default=True,
)
@click.argument('fa_in', type=click.Path(exists=True))
@click.argument('fa_out', type=click.Path(exists=False))
def filter_het(fa_in, fa_out, min_identity, min_overlap):
    """
    Filters redundant (homozygous) contigs out from an in put fasta file
    and generate a "homozygous" output fasta file.
    """
    a = mp.Aligner(fa_in)
    drop_seqs = set()
    # Align each sequence against the rest of contigs
    for query_name in a.seq_names:
        # Check if sequence has been added to drop list already
        if query_name not in drop_seqs:
            # Take best hit
            hits = a.map(a.seq(query_name))
            hit = next(hits)
            # Do not self-align
            if query_name == hit.ctg:
                try:
                    hit = next(hits)
                except StopIteration:
                    continue
            # If sequence identity is good enough, we have heterozygous pair
            identity = hit.mlen / hit.blen
            if identity >= min_identity:
                # Remove smallest contig of the pair if overlap is sufficient
                query_len = len(a.seq(query_name))
                if hit.ctg_len > query_len:
                    overlap = hit.blen / hit.ctg_len
                    if overlap >= min_overlap:
                        drop_seqs.add(query_name)
                        sys.stderr.write(
                            f"Dropped {query_name} ({query_len} bp): matched {hit.ctg}"
                            f" ({hit.ctg_len} bp) with {100 * round(identity)}%"
                            f" identity and {100 * round(overlap)}% overlap.\n"
                        )
                else:
                    overlap = hit.blen / len(a.seq(query_name))
                    if overlap >= min_overlap:
                        sys.stderr.write(
                            f"Dropped {hit.ctg} ({hit.ctg_len} bp): matched {query_name}"
                            f" ({query_len} bp) with {100 * round(identity)}%"
                            f" identity and {100 * round(overlap)}% overlap.\n"
                        )
                        drop_seqs.add(hit.ctg)

    # Write homozygous genome by omitting dropped contigs
    het_genome = SeqIO.parse(fa_in, "fasta")
    good_seqs_iter = (rec for rec in het_genome if rec.id not in drop_seqs)
    SeqIO.write(good_seqs_iter, fa_out, "fasta")


if __name__ == "__main__":
    filter_het()
