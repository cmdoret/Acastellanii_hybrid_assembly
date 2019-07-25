
# Use pilon for short reads polishing
rule pilon_polishing:
    input:
      r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq"),
      r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq"),
      assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
    output: join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
    params:
      bt2_index = temp(join(TMP, "01_Ac_{strain}_flye")),
      alignment = temp(join(TMP, 'alignments', "01_Ac_{strain}_flye.sam")),
      bt2_preset = config['params']['bowtie2'],
      pilon_preset = config['params']['pilon'],
      pilon_outdir = directory(join(TMP, 'pilon', '01_Ac_{strain}_flye'))
    singularity: "docker://quay.io/biocontainers/pilon"
    threads: 12
    shell:
      """
      bowtie2-build {input.assembly} {params.bt2_index}
      bowtie2 -x {params.bt2_index} \
              -1 {input.r1} \
              -2 {input.r2} \
              -S {params.alignment}
              -p {threads} \
              {params.bt2_preset}

      pilon --frags {params.alignment} \
            --genome {input.assembly} \
            --outdir {params.pilon_outdir} \
            --output {wildcards.strain} \
            {params.pilon_preset}

      cp {params.pilon_outdir}/
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
    singularity: "docker://quay.io/biocontainers/racon"
    shell: 
      """
      racon_wrapper --subsample 50000000 50 \
                    -u -t {threads} \
                    {input.illumina} \
                    c3_illumina_vs_03_pilon.sam \
                    03_Ac_c3_pilon.fa \
                    > {output}
      """
