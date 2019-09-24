
# Use flye for de-novo long reads assembly
rule flye_assembly:
  input: join(TMP, 'reads', '{strain}_long_reads_corrected.fa'),
  output: assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  params:
    flye_dir = directory(join(TMP, '{strain}', 'flye'))
  threads: CPUS
  singularity: "docker://cmdoret/flye:2.3.6"
  shell:
    """
    mkdir -p {params.flye_dir}
    flye --nano-corr {input} \
         --threads {threads} \
         --iterations 3 \
         -o {params.flye_dir} \
         -g 45m
    mv {params.flye_dir}/scaffolds.fasta {output}
    """


