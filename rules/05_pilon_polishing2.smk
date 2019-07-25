
# Use 2 rounds of pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end1.fq"),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  output: join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa')
  params:
    alignment = temp(join(TMP, "alignments", "05_Ac_{strain}_instagraal.sam")),
    bt2_index = temp(join(TMP, "05_Ac_{strain}_instagraal")),
    bt2_preset = config['params']['bowtie2'],
    pilon_preset = config['params']['pilon'],
    pilon_outdir = temp(join(TMP, "pilon", "05_Ac_{strain}_instagraal")),
  singularity: "docker://quay.io/biocontainers/pilon"
  threads: 12
  shell:
    """
    bowtie2-build {input.assembly} {params.bt2_index}
    bowtie2 -x {params.bt2_index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -S {params.alignment} \
            -p {threads} \
            {params.bt2_preset}

    pilon --genome {input.assembly} \
          --frags {params.alignment} \
          --threads {threads} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain} \
          {params.pilon_preset}
    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """

