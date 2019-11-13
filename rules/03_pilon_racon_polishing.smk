
# Align shotgun reads to first assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '02_Ac_{strain}_consent.fa')
  output: temporary(join(TMP, 'alignments', "02_Ac_{strain}_consent.sam"))
  params:
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"),
    bt2_index = temporary(join(TMP, "02_Ac_{strain}_consent")),
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
  input: join(TMP, 'alignments', '02_Ac_{strain}_consent.sam')
  output:
    bam=temporary(join(TMP, 'alignments', '02_Ac_{strain}_consent.bam')),
    bai=temporary(join(TMP, 'alignments', '02_Ac_{strain}_consent.bam.bai'))
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
    bam = join(TMP, "alignments", "02_Ac_{strain}_consent.bam"),
    bai = join(TMP, "alignments", "02_Ac_{strain}_consent.bam.bai"),
    alignment = join(TMP, "alignments", "02_Ac_{strain}_consent.bam"),
    assembly = join(OUT, 'assemblies', '02_Ac_{strain}_consent.fa')
  output: join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
  log: join('logs', '03_pilon_polishing_{strain}.log')
  params:
    pilon_preset = config['params']['pilon'],
    pilon_outdir = directory(join(TMP, 'pilon', '02_Ac_{strain}_consent'))
  singularity: "docker://cmdoret/pilon:1.22"
  threads: CPUS
  shell:
    """
    pilon --frags {input.bam} \
          --genome {input.assembly} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain} \
          {params.pilon_preset} \
          2> {log}

    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
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
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
  output: temporary(join(TMP, "alignments", "03_Ac_{strain}_pilon_merged_shotgun.sam"))
  params:
    bt2_index = join(TMP, "03_Ac_{strain}_pilon"),
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
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa'),
    illumina = join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"),
    alignment = join(TMP, "alignments", "03_Ac_{strain}_pilon_merged_shotgun.sam")
  output: join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')
  log: join('logs', '03_racon_polishing_{strain}.log')
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
                  > {output} \
                  2> {log}
    """
