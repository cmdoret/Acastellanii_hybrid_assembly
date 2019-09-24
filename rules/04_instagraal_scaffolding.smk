# Use redundans reduction step to remove heterozygous contigs before scaffolding
rule filter_het_contigs:
  input:
    genome = join(OUT, 'assemblies', '03_Ac_{strain}_racon.fa'),
  output:
    fa = join(OUT, 'assemblies', '04_Ac_{strain}_redundans.fa'), 
    redundans = join(TMP, '{strain}_drop_haplotypes.fa')
  params:
    similarity = 90
  singularity: "docker://lpryszcz/redundans"
  shell:
    """
    /redundans/redundans.py \
      --noscaffolding \
      --nocleaning \
      --norearrangements \
      --nogapclosing \
      -v -f test/contigs.fa -o {output.redundans}"
    """

# Generate Hi-C matrix from raw reads
rule hicstuff_hic_processing:
    input:
      r1 = join(TMP, "reads", "{strain}_hic.end1.fq.gz"),
      assembly = join(OUT, 'assemblies', '04_Ac_{strain}_redundans.fa')
    output: directory(join(TMP, "hicstuff", "{strain}"))
    threads: CPUS
    resources: mem=32000
    params:
      idx = temporary(join(TMP, '04_Ac_{strain}_redundans')),
      enzyme = "DpnII",
      r2 = join(TMP, "reads", "{strain}_hic.end2.fq.gz"),
    singularity: "docker://koszullab/hicstuff:latest"
    shell:
      """
      bowtie2-build {input.assembly} {params.idx}
      hicstuff pipeline -t {threads} \
                        -e {params.enzyme} \
                        -g {params.idx} \
                        {input.r1} {params.r2} \
                        -o {output}
      """


# Perform Hi-C based scaffolding using instagraal
rule instagraal_scaffolding:
  input:
    assembly = join(OUT, 'assemblies', '04_Ac_{strain}_redundans.fa'),
    hicstuff_dir = join(TMP, "hicstuff", "{strain}")
  output: join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  params:
    instagraal_input_dir = join(TMP, 'instagraal', '{strain}')
  shell:
    """
    mkdir -p {params.instagraal_input_dir}
    cp "{input.hicstuff_dir}/abs_fragments_contacts_weighted.txt" \
       "{input.hicstuff_dir}/fragments_list.txt" \
       "{input.hicstuff_dir}/info_contigs.txt" \
       "{params.instagraal_input_dir}"
    instagraal {params.instagraal_input_dir} {input.assembly} {output}
    """
