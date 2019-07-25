
# Use 2 rounds of pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing:
  input:
    unpack(lambda wildcards: get_fastqs(wildcards, libtype='shotgun')),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  output: join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa')
  params:
    alignment = temp(join(TMP, "alignments", "05_Ac_{strain}_instagraal.sam")),
    bt2_index = temp(join(TMP, "05_Ac_{strain}_instagraal")),
    bt2_preset = config['params']['bowtie2'],
    pilon_preset = config['params']['pilon']
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
          --output {strain} \
          {params.pilon_preset}
    mv {params.pilon_outdir}/{strain}.fasta {output}
    """

