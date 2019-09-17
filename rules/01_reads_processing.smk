

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

# Split reads into two files to allow correction with reasonable resources
N_SPLITS = 4
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)]

rule split_fasta:
  input: join(TMP, 'reads', '{strain}_long_reads.fa')
  output: expand(join(TMP, 'split_reads', '{{strain}}_long_reads', '{{strain}}_long_reads.{split}.fa'), split=split_names)
  message: "Split {input} into: {output}"
  params:
    n_splits = N_SPLITS,
    split_dir = lambda w: join(TMP, 'split_reads', f'{w.strain}_long_reads')
  shell:
    """
    mkdir -p {params.split_dir}
    seqkit split2 -p {params.n_splits} \
                  -w 0 \
                  -f \
                  -1 {input} \
                  -O {params.split_dir}
    """


# Use CONSENT to correct ONT reads
rule long_reads_correction:
  input: join(TMP, 'split_reads', '{strain}_long_reads', '{strain}_long_reads.{split}.fa')
  output: join(TMP, 'split_reads', '{strain}_long_reads_corrected.{split}.fa')
  message: "Using CONSENT to orrect {input} into {output}"
  params:
    tmp = join(TMP, "CONSENT")
  threads: 54
  singularity: "docker://cmdoret/consent:latest"
  shell:
    """
    CONSENT-correct --in {input} \
                    --tmpdir {params.tmp} \
                    --out {output} \
                    --type ONT \
                    --nproc {threads}
    """

rule merge_corrected:
  input: expand(join(TMP, 'split_reads', '{{strain}}_long_reads_corrected.{split}.fa'), split=split_names)
  output: join(TMP, 'reads', '{strain}_long_reads_corrected.fa')
  message: "Merge CONSENT outputs: {input} into {output}"
  shell: "cat {input} > {output}"
