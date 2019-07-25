
# Use pilon for short reads polishing
rule pilon_polishing:
    input:
      unpack(lambda wildcards: get_fastqs(wildcards, libtype='shotgun')),
      assembly = join(OUT, 'assemblies', '02_Ac_{strain}_haplofilter.fa')
    output: join(OUT, 'assemblies', '03_Ac_{strain}_pilon.fa')
    params:
      bt2_index = temp(join(TMP, "02_Ac_{strain}_haplofilter")),
      alignment = temp(join(TMP, 'alignments', "02_Ac_{strain}.sam")),
      bt2_preset = config['params']['bowtie2'],
      pilon_preset = config['params']['pilon'],
      pilon_outdir = directory(join(TMP, 'pilon', '02_{strain}_haplofilter'))
    singularity: "docker://quay.io/biocontainers/pilon"
    threads: 12
    shell:
      """
      bowtie2-build {input.assembly} {params.bt2_index}
      bowtie2 -x {params.bt2_index} \
              -1 {input.shotgun1} \
              -2 {input.shotgun2} \
              -S {params.alignment}
              -p {threads} \
              {params.bt2_preset}

      pilon --frags {params.alignment} \
            --genome {input.assembly} \
            --outdir {params.pilon_outdir} \
            --output {strain} \
            {params.pilon_preset}

      cp {params.pilon_outdir}/
      """

rule combine_fq_pairs:
    input: unpack(lambda wildcards: get_fastqs(wildcards, libtype='shotgun'))
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
                    PM65_merged.fq.gz \
                    c3_illumina_vs_03_pilon.sam \
                    03_Ac_c3_pilon.fa \
                    > Acastellanii_c3_pilon_racon.fasta
      """
