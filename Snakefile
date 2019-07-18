# Docker-based pipeline for Acanthamoeba castellanii genome assembly orchestrated using snakemake
# cmdoret, 20190718

strains = ["Neff", "C3"]

rule all:
  input: expand(join('out', '{strain}', 'fa', 'Ac_{strain}_pilon,fa'), strain=strains)

# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join('in', '{strain}', 'fq', 'long_reads', '{strain}_ONT_merged.fq')
  output: join('tmp', '{strain}', '{strain}_ONT_polished.fq'

# Use flye for de-novo long reads assembly
rule flye_assembly:
  input:
    join('tmp', '{strain}', '{strain}_ONT_polished.fq'),
  output: join('out', '{strain}', 'fa', '01_Ac_{strain}_flye.fa')

# Use mummer to compute all-vs-all contig similarity and filter out small contigs 
# of unmerged haplotypes
rule filter_haplotypes:
  input: join('out', '{strain}', 'fa', '01_Ac_{strain}_flye.fa')
  output:
    keep = join('out', '{strain}', 'fa', '02_Ac_{strain}_haplofilter.fa'), 
    drop = join('tmp', '{strain}', 'drop_haplotypes.fa')
  params:
    similarity = 90

# Use pilon for short reads polishing
rule pilon_polishing:
  input:
    assembly = join('out', '{strain}', 'fa', '02_Ac_{strain}_haplofilter.fa'),
    illumina = join('in', '{strain}', 'fq', '{strain}_shotgun.fq')
  output: join('out', '{strain}', 'fa', '03_Ac_{strain}_pilon.fa')

# Use racon-illumina for another round of short reads polishing
rule racon_polishing:
  input:
    assembly = join('out', '{strain}', 'fa', '03_Ac_{strain}_pilon.fa'),
    illumina = join()
  output: join('out', '{strain}', 'fa', '04_Ac_{strain}_racon.fa'

# Generate Hi-C matrix from raw reads
rule: hicstuff_hic_processing:
  input:
    assembly = join(),
    hic_end1 = join(),
    hic_end2 = join()
  output: join()

# Perform Hi-C based scaffolding using instagraal
rule instagraal_scaffolding:
  input:
    assembly = join(),
    matrix = join()
  output: join('out', '{strain}', 'fa', '05_Ac_{strain}_instagraal.fa')

rule post_hic_pilon_polishing:
  input:
    assembly = join(),
    illumina = join()
  output: join('out', '{strain}', 'fa', '06_Ac_{strain}_pilon.fa')
