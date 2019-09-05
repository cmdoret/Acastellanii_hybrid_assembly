
# Align shotgun reads to first assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = GS.remote(join(TMP, "reads", "{strain}_shotgun.end1.fq.gz")),
    r2 = GS.remote(join(TMP, "reads", "{strain}_shotgun.end2.fq.gz")),
    assembly = GS.remote(join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa'))
  output: GS.remote(join(TMP, 'alignments', "01_Ac_{strain}_flye.sam"))
  params:
    bt2_index = temp(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: CPUS
  resources: mem=32000
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {output} \
            -p {threads} \
            {params.bt2_preset}
    """

# Use pilon for short reads polishing of the long reads assembly
rule pilon_polishing:
  input:
    alignment = GS.remote(join(TMP, "alignments", "01_Ac_{strain}_flye.sam")),
    assembly = GS.remote(join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa'))
  output: GS.remote(join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa'))
  params:
    pilon_preset = config['params']['pilon'],
    pilon_outdir = GS.remote(directory(join(TMP, 'pilon', '01_Ac_{strain}_flye')))
  singularity: "docker://cmdoret/pilon:1.22"
  threads: CPUS
  resources: mem=256000
  shell:
    """
    pilon --frags {input.alignment} \
          --genome {input.assembly} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain} \
          {params.pilon_preset}

    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """

rule combine_fq_pairs:
  input:
    r1 = GS.remote(join(TMP, "reads", "{strain}_shotgun.end1.fq.gz")),
    r2 = GS.remote(join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"))
  output: GS.remote(join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"))
  params:
  shell: "cat {input.r1} {input.r2} > {output}"

# Map shotgun reads in single end mode for consensus correction via racon
rule align_merged_shotgun_pilon_assembly:
  input:
    reads = GS.remote(join(TMP, "reads", "{strain}_merged_shotgun.fq.gz")),
    assembly = GS.remote(join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa'))
  output: GS.remote(join(TMP, "alignments", "03_Ac_{strain}_pilon_merged_shotgun.sam"))
  params:
    bt2_index = GS.remote(join(TMP, "01_Ac_{strain}_flye")),
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
    assembly = GS.remote(join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')),
    illumina = GS.remote(join(TMP, "reads", "{strain}_merged_shotgun.fq.gz")),
    alignment = GS.remote(join(TMP, "alignments", "03_Ac_{strain}_pilon_merged_shotgun.sam"))
  output: GS.remote(join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa'))
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
