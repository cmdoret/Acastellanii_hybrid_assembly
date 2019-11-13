
# Use flye for de-novo long reads assembly
rule flye_assembly:
  input: join(TMP, 'reads', '{strain}_long_reads_filtered.fa'),
  output: assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  params:
    flye_dir = directory(join(TMP, '{strain}', 'flye'))
  threads: CPUS
  singularity: "docker://cmdoret/flye:2.3.6"
  shell:
    """
    mkdir -p {params.flye_dir}
    flye --nano-raw {input} \
    --threads {threads} \
    --iterations 3 \
    -o {params.flye_dir} \
    -g 45m
    mv {params.flye_dir}/scaffolds.fasta {output}
    """

# Polish the draft with long reads
rule long_reads_polishing:
  input:
    ont = join(TMP, 'reads', '{strain}_long_reads.fa'),
    contigs = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: join(OUT, 'assemblies', '02_Ac_{strain}_consent.fa')
  message: "Using CONSENT to polish the draft assembly with ONT reads."
  params:
    tmp = join(TMP, "CONSENT")
  threads: 54
  singularity: "docker://cmdoret/consent:1.2"
  shell:
    """
    CONSENT-polish  --contigs {input.contigs} \
                    --reads {input.ont} \
                    --out {output} \
                    --nproc {threads}
    """
