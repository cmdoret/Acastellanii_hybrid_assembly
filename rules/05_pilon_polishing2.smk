
#Align shotgun reads to the scaffolded assembly
rule align_shotgun_3c_assembly:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish.fa')
  output: temporary(join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.sam"))
  params:
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"),
    bt2_index = temporary(join(TMP, "04_Ac_{strain}_instagraal_polish")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  conda: '../envs/align.yaml'
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

rule index_shotgun_bam_pilon_3c:
  input: join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.sam")
  output: 
    bam = temporary(join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.bam")),
    bai = temporary(join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.bam.bai"))
  singularity: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
  conda: '../envs/align.yaml'
  threads: CPUS
  shell:
    """
    samtools sort -@ {threads} -O BAM -o {output.bam} {input} 
    samtools index -@ {threads} {output.bam}
    """

# Use pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing:
  input:
    bai = join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.bam.bai"),
    alignment = join(TMP, "alignments", "04_Ac_{strain}_instagraal_polish.bam"),
    assembly = join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish.fa')
  output: join(OUT, 'assemblies', '05_Ac_{strain}_pilon.fa')
  params:
    #pilon_preset = config['params']['pilon'],
    pilon_outdir = temp(join(TMP, "pilon", "04_Ac_{strain}_instagraal_polish"))
  singularity: "quay.io/biocontainers/pilon:1.22--py36_0"
  conda: '../envs/pilon.yaml'
  threads: CPUS
  resources: mem=256000
  shell:
    """
    pilon --genome {input.assembly} \
          --frags {input.alignment} \
          --threads {threads} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain}
    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """

# Align shotgun reads to the polished assembly
rule align_shotgun_3c_assembly_round2:
  input:
    r1 = join(TMP, "reads", "{strain}_shotgun.end1.fq.gz"),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_pilon.fa')
  output: temp(join(TMP, "alignments", "05_Ac_{strain}_pilon.sam"))
  params:
    r2 = join(TMP, "reads", "{strain}_shotgun.end2.fq.gz"),
    bt2_index = temp(join(TMP, "05_Ac_{strain}_pilon")),
    bt2_preset = config['params']['bowtie2']
  singularity: "docker://cmdoret/bowtie2:2.3.4.1"
  conda: '../envs/align.yaml'
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

rule index_shotgun_bam_pilon_3c_round2:
  input: join(TMP, "alignments", "05_Ac_{strain}_pilon.sam")
  output: 
    bam = temporary(join(TMP, "alignments", "05_Ac_{strain}_pilon.bam")),
    bai = temporary(join(TMP, "alignments", "05_Ac_{strain}_pilon.bam.bai"))
  singularity: "docker://biocontainers/samtools:v1.7.0_cv4"
  conda: '../envs/align.yaml'
  threads: CPUS
  shell:
    """
    samtools sort -@ {threads} -O BAM -o {output.bam} {input} 
    samtools index -@ {threads} {output.bam}
    """

# Use pilon to polish the HI-C scaffolded assembly
rule post_hic_pilon_polishing_round2:
  input:
    alignment = join(TMP, "alignments", "05_Ac_{strain}_pilon.bam"),
    assembly = join(OUT, 'assemblies', '05_Ac_{strain}_pilon.fa'),
    bai = join(TMP, "alignments", "05_Ac_{strain}_pilon.bam.bai")
  output: join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa')
  params:
    #pilon_preset = config['params']['pilon'],
    pilon_outdir = temp(join(TMP, "pilon", "05_Ac_{strain}_pilon"))
  singularity: "quay.io/biocontainers/pilon:1.22--py36_0"
  conda: '../envs/pilon.yaml'
  threads: CPUS
  resources: mem=256000
  shell:
    """
    pilon --genome {input.assembly} \
          --frags {input.alignment} \
          --threads {threads} \
          --outdir {params.pilon_outdir} \
          --output {wildcards.strain}
    mv {params.pilon_outdir}/{wildcards.strain}.fasta {output}
    """

