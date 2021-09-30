from typing import List

def get_fq_units(wildcards, end: int=1) -> List[str]:
    """
    Get fastq files (units) of a particular library type of one sample
    from the unit sheet. end can be 1 or 2 (i.e. forward or reverse).
    """
    fqs = units.loc[
      (units.strain == wildcards.strain) &
      (units.libtype == wildcards.libtype),
      f"fq{end}"
    ].dropna()
    return list(fqs)


# Combine all fastq files from the same sample / library type combination
# If there is a single fastq, symlink it to spare memory
rule combine_units:
  input: lambda w: get_fq_units(w, 1)
  params:
    r2 = lambda w: get_fq_units(w, 2)
  message:
    """
    Merging :
      {input} into:  {output.r1}
      {params.r2} into:  {output.r2}
    """
  output:
    r1 = join(TMP, "reads", "{strain}_{libtype}.end1.fq.gz"),
    r2 = join(TMP, "reads", "{strain}_{libtype}.end2.fq.gz")
  threads: 1
  run:
    if len(input[:]) > 1:
      shell(f"cat {' '.join(input[:])} > {output['r1']}")
      if len(params['r2']):
        shell(f"cat {' '.join(params['r2'])} > {output['r2']}")
    else:
      shell(f"ln -s {input[0]} {output['r1']}")
      if len(params['r2']):
        shell(f"ln -s {params['r2'][0]} {output['r2']}")



# convert fastq reads to fasta for correction
rule fastq_to_fasta_ONT:
  input: 
    r1 = join(TMP, "reads", "{strain}_long_reads.end1.fq.gz")
  output: join(TMP, 'reads', '{strain}_long_reads.fa')
  singularity: "docker://cmdoret/seqtk:1.3"
  conda: '../envs/seqtk.yaml'
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
    keep=60
  conda: "../envs/filtlong.yaml"
  singularity: "docker://nanozoo/filtlong:0.2.0--0c4cbe3"
  shell: "filtlong -p {params.keep} -1 {input.ilm1} -2 {input.ilm2} {input.ont} > {output}"

# Format long reads fasta into 1 read / line for CONSENT
rule format_long_reads:
  input: join(TMP, 'reads', '{strain}_long_reads_filtered_tmp.fa')
  output: join(TMP, 'reads', '{strain}_long_reads_filtered.fa')
  shell: "sed 's/^\\(>[^ ]*\\) .*/\\1/' {input} > {output}"
