

# Use ra for de-novo long reads assembly
rule ra_assembly:
  input: join(TMP, 'reads', '{strain}_long_reads_filtered.fa')
  output: join(OUT, 'assemblies', '01_Ac_{strain}_ra.fa')
  threads: 12
  singularity: 'docker://cmdoret/ra:0.2.1'
  shell: "ra -t {threads} {input} > {output}"

# Polish the draft with long reads
rule long_reads_polishing:
  input:
    ont = join(TMP, 'reads', '{strain}_long_reads.fa'),
    contigs = join(OUT, 'assemblies', '01_Ac_{strain}_ra.fa')
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
