# Use flye for de-novo long reads assembly
rule flye_assembly:
  input: join(TMP, 'reads', '{strain}_long_reads_filtered.fa')
  output: join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  log: join('logs', '02_flye_assembly_{strain}.log')
  params:
    flye_dir = directory(join(TMP, '{strain}', 'flye'))
  threads: CPUS
  singularity: "quay.io/biocontainers/flye:2.8.2"
  conda: "../envs/flye.yaml"
  shell:
    """
    mkdir -p {params.flye_dir}
    flye --nano-raw {input} \
    --threads {threads} \
    --iterations 3 \
    -o {params.flye_dir} \
    -g 45m \
    2> {log}
    
    mv {params.flye_dir}/assembly.fasta {output}
    """
