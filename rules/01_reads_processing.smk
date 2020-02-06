

rule combine_units:
  input: unpack(get_fastqs)
  output:
    r1 = join(TMP, "reads", "{strain}_{libtype}.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_{libtype}.end2.fq.gz")
  run:
    shell("cat {i1} > {o1}".format(i1=input['r1'], o1=output['r1']))
    if len(input.keys()) == 2:
      shell("cat {i2} > {o2}".format(i2=input['r2'], o2=output['r2']))
    else:
        shell("touch {o2}".format(o2=output['r2']))

# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: 
    r1 = join(TMP, "reads", "{strain}_long_reads.end1.fq.gz")
  output: join(TMP, 'reads', '{strain}_long_reads.fa')
  singularity: "docker://cmdoret/seqtk:1.3"
  shell:
    """
    seqtk seq -a {input.r1} > {output}
    """


rule filter_long_reads:
  input: 
    ont = join(TMP, 'reads', '{strain}_long_reads.fa'),
    ilm1 = join(TMP, 'reads', '{strain}_shotgun.end1.fq.gz'),
    ilm2 = join(TMP, 'reads', '{strain}_shotgun.end2.fq.gz')
  output: temporary(join(TMP, 'reads', '{strain}_long_reads_filtered_tmp.fa'))
  params:
    keep=80
  conda: "../envs/filtlong.yaml"
  shell: "filtlong -p {params.keep} -1 {input.ilm1} -2 {input.ilm2} {input.ont} > {output}"

# Format long reads fasta into 1 read / line for CONSENT
rule format_long_reads:
  input: join(TMP, 'reads', '{strain}_long_reads_filtered_tmp.fa')
  output: join(TMP, 'reads', '{strain}_long_reads_filtered.fa')
  shell: "sed 's/^\\(>[^ ]*\\) .*/\\1/' {input} > {output}"
