
#Align shotgun reads to the scaffolded assembly
rule align_shotgun_3c_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  output: temp(join(TMP, "alignments", "05_Ac_{strain}_instagraal.sam"))
  params:
    bt2_index = temp(join(TMP, "05_Ac_{strain}_instagraal")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: 12
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

# Use pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing:
  input:
    alignment = join(TMP, "alignments", "05_Ac_{strain}_instagraal.sam"),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  output: join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa')
  params:
    pilon_preset = config['params']['pilon'],
    pilon_outdir = temp(join(TMP, "pilon", "05_Ac_{strain}_instagraal")),
  singularity: "docker://cmdoret/pilon:1.22"
  threads: 12
  resources: mem=256000
  shell:
    """
    pilon --genome {input.assembly} \
          --frags {input.alignment} \
          --threads {threads} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain} \
          {params.pilon_preset}
    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """

