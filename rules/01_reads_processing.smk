

# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: unpack(lambda wildcards: get_fastqs(wildcards, libtype='ONT'))
  output: join(TMP, 'reads', '{strain}_ONT.fa')
  singularity: "docker://alpine"
  shell:
    """
    paste - - - - \
      < {input.r1} \
      | cut -f 1,2 \
      | sed 's/^@/>/' \
      | tr "\t" "\n" \
      > {output}
    """


# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join(TMP, 'reads', '{strain}_ONT.fa')
  output: join(TMP, 'reads', '{strain}_ONT_polished.fa')
  singularity: "docker://cmdoret/consent"
  shell: "CONSENT-correct --in {input} --out {output} --type ONT"

