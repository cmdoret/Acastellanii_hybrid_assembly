
rule combine_units:
  input:
    lambda w: GS.remote(units_dict[f"{w.libtype}"][f"{w.strain}"][f"fq{w.end}"]),
  message:
    """
    Merging :
      {input} into:  {output}
    """
  output:
    GS.remote(join(TMP, "reads", "{strain}_{libtype}.end{end}.fq.gz")),
  threads: 12
  shell: "cat {input} > {output}"


# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: 
    r1 = GS.remote(join(TMP, "reads", "{strain}_long_reads.end1.fq.gz"))
  output: GS.remote(join(TMP, 'reads', '{strain}_long_reads.fa'))
  singularity: "docker://cmdoret/seqtk:1.3"
  shell:
    """
    seqtk seq -a {input.r1} > {output}
    """


# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: GS.remote(join(TMP, 'reads', '{strain}_long_reads.fa'))
  output: GS.remote(join(TMP, 'reads', '{strain}_long_reads_corrected.fa'))
  params:
    tmp = join(TMP, "CONSENT")
  threads: 5
  singularity: "docker://cmdoret/consent:latest"
  shell:
    """
    /app/CONSENT-correct --in {input} \
                    --tmpdir {params.tmp} \
                    --out {output} \
                    --type ONT \
                    --nproc {threads}
    """
