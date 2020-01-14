
# Align shotgun reads to the ONT assembly
rule align_shotgun_ont_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: temporary(join(TMP, 'alignments', "shotgun_vs_01_Ac_{strain}_flye.bam"))
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
            -p {threads} \
            {params.bt2_preset} \
      | samtools sort -@ {threads} -o {output}
    """


# Align ONT reads to back to the ONT assembly
rule align_ont_ont_assembly:
  input:
    ont = join(TMP, "reads", "{strain}_long_reads.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: temporary(join(TMP, 'alignments', "ONT_vs_01_Ac_{strain}_flye.bam"))
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  threads: CPUS
  resources: mem=32000
  shell:
    """
    minimap2 -ax map-ont \
            {input.assembly} {input.ont} \
            -t {threads} \
      | samtools sort -@ {threads} -o {output}
    """

# Use pilon for short reads polishing of the long reads assembly
# TODO: Estimate genome size and coverage from input
rule hypo_polishing:
  input:
    bam_sr = join(TMP, "alignments", "shotgun_vs_01_Ac_{strain}_flye.bam"),
    bam_ont = join(TMP, "alignments", "ONT_vs_01_Ac_{strain}_flye.bam"),
    assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output:
    assembly = join(OUT, 'assemblies', '02_Ac_{strain}_hypo.fa'),
    fqlist = temporary(join(TMP, '{strain}_fq_list.txt'))
  log: join('logs', '02_hypo_polishing_{strain}.log')
  singularity: "docker://cmdoret/pilon:1.22"
  params:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz")
  threads: CPUS
  shell:
    """
    echo -e "{params.r1}\n{params.r2}" > {output.fqlist}
    hypo -t {threads} \
         -o {output.assembly} \
         -r "@"{output.fqlist} \
         -b {input.bam_sr} \
         -B {input.bam_ont} \
         -d {input.assembly} \
         -s 45m \
         -c 60 \
          2> {log}
    """
