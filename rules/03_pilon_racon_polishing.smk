
# Align shotgun reads to first assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: temp(join(TMP, 'alignments', "01_Ac_{strain}_flye.sam"))
  params:
    bt2_index = temp(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: 12
  resources: mem="32G"
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {output.alignment}
            -p {threads} \
            {params.bt2_preset}
    """

# Use pilon for short reads polishing of the long reads assembly
rule pilon_polishing:
  input:
    alignment = join(TMP, "alignments", "01_Ac_{strain}_flye.sam"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
  params:
    pilon_preset = config['params']['pilon'],
    pilon_outdir = directory(join(TMP, 'pilon', '01_Ac_{strain}_flye'))
  singularity: "docker://cmdoret/pilon:1.22"
  threads: 12
  resources: mem="256G"
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
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz")
  output: temp(join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"))
  shell: "cat {input.r1} {input.r2} > {output}"

# Map shotgun reads in single end mode for consensus correction via racon
rule align_merged_shotgun_pilon_assembly:
  input:
    reads = join(TMP, "reads", "{strain}_merged_shotgun.fq.gz"),
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
  output: temp(join(TMP, "alignments", "03_Ac_{strain}_pilon_merged_shotgun.sam"))
  params:
    bt2_index = temp(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: 12
  resources: mem="32G"
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -U {input.reads} \
            -S {output.alignment}
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
  threads: 12
  resources: mem="128G"
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
