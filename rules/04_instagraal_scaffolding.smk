
# Generate Hi-C matrix of the nuclear genome from raw reads 
rule hicstuff_hic_processing:
    input:
      r1 = join(TMP, "reads", "{strain}_hic.end1.fq.gz"),
      assembly = join(OUT, 'assemblies', '02_Ac_{strain}_hypo_nucl.fa')
    output: directory(join(TMP, "hicstuff", "{strain}_nucl"))
    threads: CPUS
    resources: mem=32000
    params:
      idx = temporary(join(TMP, '02_Ac_{strain}_hypo_nucl')),
      enzyme = "DpnII",
      r2 = join(TMP, "reads", "{strain}_hic.end2.fq.gz"),
    singularity: "docker://koszullab/hicstuff:v2.2.2"
    shell:
      """
      bowtie2-build {input.assembly} {params.idx}
      hicstuff pipeline -t {threads} \
                        -e {params.enzyme} \
                        -g {params.idx} \
                        {input.r1} {params.r2} \
                        -o {output} \
                        -T {output}
      """


# Perform Hi-C based scaffolding using instagraal
rule instagraal_scaffolding:
  input:
    assembly = join(OUT, 'assemblies', '02_Ac_{strain}_hypo_nucl.fa'),
    hicstuff_dir = join(TMP, "hicstuff", "{strain}_nucl")
  output:
    assembly = join(OUT, 'assemblies', '03_Ac_{strain}_instagraal_nucl.fa'),
    frags = join(TMP, 'info_frags_{strain}.txt')
  params:
    instagraal_outdir = join(TMP, 'instagraal', '{strain}_nucl')
  shell:
    """
    # Run instagraal
    instagraal {input.hicstuff_dir} {input.assembly} {params.instagraal_outdir}
    # Take the assembly out of instagraal output dir
    cp {params.instagraal_outdir}/$(basename {input.hicstuff_dir})/test_mcmc_4/genome.fasta {output.assembly}
    # Take info_frags file out, it will be useful for instagraal-polish
    cp {params.instagraal_outdir}/$(basename {input.hicstuff_dir})/test_mcmc_4/info_frags.txt {output.frags}
    """


rule instagraal_polish:
  input:
    fasta = join(OUT, 'assemblies', '02_Ac_{strain}_hypo_nucl.fa'),
    frags = join(TMP, 'info_frags_{strain}.txt')
  output: join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish_nucl.fa')
  shell: "instagraal-polish -m polishing -i {input.frags} -f {input.fasta} -o {output}"


# Merge nuclear scaffolds and mitochondrial contigs
rule merge_organelles:
    input:
        nucl = join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish_nucl.fa'),
        mito = join(OUT, 'assemblies', '02_Ac_{strain}_hypo_mito.fa')
    output: join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish.fa')
    shell: "cat {input.nucl} {input.mito} > {output}"
