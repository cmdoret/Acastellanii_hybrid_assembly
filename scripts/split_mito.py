# This script uses minimap2 python binding to separate mitochondrial from nuclear contigs
# in the input assembly. It does so by aligning contigs to the mitochondrial genome from
# a previous assembly and applying a threshold for alignment identity and proportion of aligned
# sequence.
# cmdoret, 20191031

import mappy as mp
import itertools as it
import click
from Bio import SeqIO
import sys




@click.command()
@click.option(
    "-I",
    "--min-identity",
    help="Minimum (blast-like) identity required to consider mitochondrial match",
    default=0.90,
    show_default=True,
)
@click.option(
    "-O",
    "--min-overlap",
    help="Minimum overlap required to consider mitochondrial contigs. Defined as: length(alignment) / length(c) where c is the smallest contig of the aligned pair.",
    default=0.3,
    show_default=True,
)
@click.option(
    "-g",
    "--genome",
    help="Input genome assembly.",
    type=click.Path(exists=True),
)
@click.option(
    "-m",
    "--mitochondrion",
    help="Input mitochondrial sequence used for aligning contigs.",
    type=click.Path(exists=True),
)
@click.option(
    "-N",
    "--nuclear-out",
    help="Output fasta with nuclear contigs.",
    type=click.Path(exists=False),
)
@click.option(
    "-M",
    "--mitochondrion-out",
    help="Output fasta with mitochondrial contigs.",
    type=click.Path(exists=False),
)
def split_mitochondrial(genome, mitochondrion, nuclear_out, mitochondrion_out, min_identity, min_overlap):
    """
    Split nuclear and mitochondrial contigs into two different fasta files.
    """
    s = mp.Aligner(mitochondrion, preset='asm5')
    q = mp.Aligner(genome, preset='asm5')
    mito_seqs = set()
    # Align each sequence against the rest of contigs
    for contig_name in q.seq_names:
        # Take best hit
        contig_seq = q.seq(contig_name)
        contig_len = len(contig_seq)
        hits = s.map(contig_seq)
        for hit in hits:
            # If sequence identity is good enough, we have a mito contig 
            identity = hit.mlen / hit.blen
            overlap = hit.blen / contig_len
            if identity >= min_identity:
                if overlap >= min_overlap:
                    mito_seqs.add(contig_name)
                    sys.stderr.write(
                        f"Detected mitochondrial sequence {contig_name} ({contig_len} bp)"
                        f" with {100 * round(identity, 2)}%"
                        f" identity and {100 * round(overlap, 2)}% overlap.\n"
                    )
    # Write output fastas
    with open(mitochondrion_out, 'w') as mito_out, open(nuclear_out, 'w') as nucl_out:
        for record in SeqIO.parse(genome, 'fasta'):
            if record.id in mito_seqs:
                SeqIO.write(record, mito_out, 'fasta')
            else:
                SeqIO.write(record, nucl_out, 'fasta')

if __name__ == "__main__":
    split_mitochondrial()
