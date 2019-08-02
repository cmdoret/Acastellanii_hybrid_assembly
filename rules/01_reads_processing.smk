

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


# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join(TMP, 'reads', '{strain}_long_reads.fa')
  output: join(TMP, 'reads', '{strain}_long_reads_corrected.fa')
  params:
    tmp = join(TMP, "CONSENT")
  threads: 12
  singularity: "docker://cmdoret/consent:latest"
  shell:
    """
    /app/CONSENT-correct --in {input} \
                    --tmpdir {params.tmp} \
                    --out {output} \
                    --type ONT \
                    --nproc {threads}
    """
