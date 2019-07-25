

rule combine_units:
  input: unpack(get_fastqs)
  output:
    r1 = join(TMP, "reads", "{strain}_{libtype}.end1.fq"),
    r2 = join(TMP, "reads", "{strain}_{libtype}.end2.fq")
  run:
    print(len(input))
    print(input[:])
    shell("cat {i1} > {o1}".format(i1=input[0], o1=output['r1']))
    if len(input) == 2:
      shell("cat {i2} > {o2}".format(i2=input[1], o2=output['r2']))
    else:
        shell("touch {o2}".format(o2=output['r2']))

# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: 
    r1 = join(TMP, "reads", "{strain}_long_reads.end1.fq")
  output: join(TMP, 'reads', '{strain}_long_reads.fa')
  singularity: "docker://quay.io/biocontainers/seqtk:1.3--ha92aebf_0"
  shell:
    """
    seqtk seq -a {input.r1} > {output}
    """


# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join(TMP, 'reads', '{strain}_long_reads.fa')
  output: join(TMP, 'reads', '{strain}_long_reads_polished.fa')
  singularity: "docker://cmdoret/consent"
  shell: "CONSENT-correct --in {input} --out {output} --type ONT"

