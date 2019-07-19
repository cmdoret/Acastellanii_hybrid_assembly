# Docker-based pipeline for Acanthamoeba castellanii genome assembly orchestrated using snakemake
# cmdoret, 20190718

from os.path import join

strains = ["Neff", "C3"]
IN = join('data', IN)
OUT = join('data', OUT)
TMP = join('data', TMP)

rule all:
  input: expand(join(OUT, '{strain}', 'assemblies', 'Ac_{strain}_pilon,fa'), strain=strains)

# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: join(IN, '{strain}', 'reads', 'long_reads', '{strain}_ONT_merged,fq')
  output: join(IN, '{strain}', 'reads', 'long_reads', '{strain}_ONT_merged,fa')
  singularity: docker://alpine
  shell:
    """
    paste - - - - \
      < {input} \
      | cut -f 1,2 \
      | sed 's/^@/>/' \
      | tr "\t" "\n" \
      > {output}
    """
# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join(IN, '{strain}', 'reads', 'long_reads', '{strain}_ONT_merged.fa')
  output: join(TMP, '{strain}', '{strain}_ONT_polished.fa'
  singularity: docker://cmdoret/consent
  shell: "CONSENT-correct --in {input} --out {output} --type ONT"

# Use flye for de-novo long reads assembly
rule flye_assembly:
  input:
    join(TMP, '{strain}', '{strain}_ONT_polished.fa'),
  output: assembly = join(OUT, '{strain}', 'assemblies', '01_Ac_{strain}_flye.fa')
  params:
    flye_dir = dir(join(TMP, '{strain}', 'flye'))
  docker pull quay.io/biocontainers/mummer:3.23--pl526_7
  threads: 12
  shell:
    """
    flye --nano-corrected {input.reads} \
         --threads {threads} \
         --iterations 3 \
         -o {params.flye_dir} \
         -g 45m
    mv {params.flye_dir}/scaffolds.fa {output}
    """


# Use mummer to compute all-vs-all contig similarity
rule dnadiff_pairwise:
  input: join(OUT, '{strain}', 'assemblies', '01_Ac_{strain}_flye.fa')
  output: join(TMP, '{strain}', 'contig_similarity.tsv')
  singularity: docker://quay.io/biocontainers/mummer
  script: "scripts/pairwise_dnadiff {input} {output}"

# filter out small contigs of unmerged haplotypes
rule filter_haplotypes:
  input: join(OUT, '{strain}', 'assemblies', '01_Ac_{strain}_flye.fa')
  output:
    keep = join(OUT, '{strain}', 'assemblies', '02_Ac_{strain}_haplofilter.fa'), 
    drop = join(TMP, '{strain}', 'drop_haplotypes.fa')
  params:
    similarity = 90
  

# Use pilon for short reads polishing
rule pilon_polishing:
  input:
    assembly = join(OUT, '{strain}', 'assemblies', '02_Ac_{strain}_haplofilter.fa'),
    illumina = join(IN, '{strain}', 'reads', '{strain}_shotgun.fq')
  output: join(OUT, '{strain}', 'fa', '03_Ac_{strain}_pilon.fa')
  singularity

# Use racon-illumina for another round of short reads polishing
rule racon_polishing:
  input:
    assembly = join(OUT, '{strain}', 'fa', '03_Ac_{strain}_pilon.fa'),
    illumina = join()
  output: join(OUT, '{strain}', 'fa', '04_Ac_{strain}_racon.fa'

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
    assembly = join(IN, '{strain}', 'fa', '04_Ac_{strain}_racon.fa'),
    matrix = dir(join(OUT, '{strain}', 'hicstuff', '04_Ac_{strain}_racon'))
  output: join(OUT, '{strain}', 'fa', '05_Ac_{strain}_instagraal.fa')

# Use 2 rounds of pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing:
  input:
    assembly = join(),
    illumina = join()
  output: join(OUT, '{strain}', 'fa', '06_Ac_{strain}_pilon2.fa')

