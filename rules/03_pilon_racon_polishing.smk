
# Align shotgun reads to first assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: temporary(join(TMP, 'alignments', "01_Ac_{strain}_flye.sam"))
  params:
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"),
    bt2_index = temporary(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: CPUS
  resources: mem=32000
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -1 {input.r1} \
            -2 {params.r2} \
            -S {output} \
            -p {threads} \
            {params.bt2_preset}
    """

rule index_shotgun_bam_pilon:
  input: join(TMP, 'alignments', '01_Ac_{strain}_flye.sam')
  output:
    bam=temporary(join(TMP, 'alignments', '01_Ac_{strain}_flye.bam')),
    bai=temporary(join(TMP, 'alignments', '01_Ac_{strain}_flye.bam.bai'))
  singularity: "docker://biocontainers/samtools:v1.7.0_cv4"
  threads: CPUS
  shell:
    """
    samtools sort -@ {threads} -O BAM -o {output.bam} {input}
    samtools index -@ {threads} {output.bam}
    """

# Use pilon for short reads polishing of the long reads assembly
rule pilon_polishing:
  input:
    bam = join(TMP, "alignments", "01_Ac_{strain}_flye.bam"),
    bai = join(TMP, "alignments", "01_Ac_{strain}_flye.bam.bai"),
    alignment = join(TMP, "alignments", "01_Ac_{strain}_flye.bam"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: join(OUT, 'assemblies', '02_Ac_{strain}_pilon.fa')
  params:
    pilon_preset = config['params']['pilon'],
    pilon_outdir = directory(join(TMP, 'pilon', '01_Ac_{strain}_flye'))
  singularity: "docker://cmdoret/pilon:1.22"
  threads: CPUS
  resources: mem=256000
  shell:
    """
    pilon --frags {input.bam} \
          --genome {input.assembly} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain} \
          {params.pilon_preset}

    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """


# Use minimap2-based custom tool to filter out redundant contigs
# Each contigs is mapped against the genome. If it has high BLAST-like
# identity and overlap with a larger contig, it is considered homozygous
# and discarded.
rule filter_het_contigs:
  input: join(OUT, 'assemblies', '02_Ac_{strain}_pilon.fa'),
  output: join(OUT, 'assemblies', '03_Ac_{strain}_homozygous.fa')
  params:
    identity = 0.51,
    overlap = 0.8
  shell:
    """
    python scripts/filter_het.py \
      -I {params.identity} \
      -O {params.overlap} \
      {input} \
      {output}
    """


rule combine_fq_pairs:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz")
  output: temporary(join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"))
  params:
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz")
  shell: "cat {input.r1} {params.r2} > {output}"


# Map shotgun reads in single end mode for consensus correction via racon
rule align_merged_shotgun_pilon_assembly:
  input:
    reads = join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"),
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_homozygous.fa')
  output: temporary(join(TMP, "alignments", "03_Ac_{strain}_homozygous_merged_shotgun.sam"))
  params:
    bt2_index = temporary(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: CPUS
  resources: mem=32000
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -U {input.reads} \
            -S {output} \
            -p {threads} \
            {params.bt2_preset}
    """

# Use racon-illumina for another round of short reads polishing
rule racon_polishing:
  input:
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_homozygous.fa'),
    illumina = join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"),
    alignment = join(TMP, "alignments", "03_Ac_{strain}_homozygous_merged_shotgun.sam")
  output: join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')
  threads: CPUS
  resources: mem=128000
  singularity: "docker://cmdoret/racon:1.3.2"
  shell: 
    """
    racon_wrapper --subsample 50000000 50 \
                  -u -t {threads} \
                  {input.illumina} \
                  {input.alignment} \
                  {input.assembly} \
                  > {output}
    """
