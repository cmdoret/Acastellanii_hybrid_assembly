
# Align shotgun reads to first assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: temp(join(TMP, 'alignments', "01_Ac_{strain}_flye.sam"))
  params:
    bt2_index = temp(join(TMP, "01_Ac_{strain}_flye")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: 12
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
  input: r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq")
  output: temp(join(TMP, "reads", "{strain}_merged_shotgun.fq"))
  shell: "cat {input} > {output}"

# Use racon-illumina for another round of short reads polishing
rule racon_polishing:
  input:
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa'),
    illumina = join(TMP, "reads", "{strain}_merged_shotgun.fq")
  output: join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')
  threads: 12
  singularity: "docker://cmdoret/racon:1.3.2"
  shell: 
    """
    racon_wrapper --subsample 50000000 50 \
                  -u -t {threads} \
                  {input.illumina} \
                  c3_illumina_vs_03_pilon.sam \
                  03_Ac_c3_pilon.fa \
                  > {output}
    """
