
# Generate Hi-C matrix from raw reads
rule hicstuff_hic_processing:
    input:
      r1 = join(TMP, "reads", "{strain}_Hi-C.end1.fq"),
      r2 = join(TMP, "reads", "{strain}_Hi-C.end1.fq"),
      assembly = join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')
    output: directory(join(TMP, "hicstuff", "{strain}"))
    threads: 12
    params:
      enzyme = "DpnII"
    singularity: "docker://koszullab/hicstuff:latest"
    shell:
      """
      hicstuff pipeline -t {threads} \
                        -e {params.enzyme} \
                        -g {input.assembly} \
                        {input.r1} {input.r2} \
                        -o {output}
      """


# Perform Hi-C based scaffolding using instagraal
rule instagraal_scaffolding:
  input:
    assembly = join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa'),
    hicstuff_dir = join(TMP, "hicstuff", "{strain}")
  output: join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa')
  params:
    instagraal_input_dir = join(TMP, 'instagraal', '{strain}')
  shell:
    """
    mkdir {params.instagraal_input_dir}
    cp "{input.hicstuff_dir}/abs_fragments_contacts_weighted.txt" \
       "{input.hicstuff_dir}/fragments_list.txt" \
       "{input.hicstuff_dir}/info_contigs.txt" \
       "{params.instagraal_input_dir}"
    instagraal {params.instagraal_input_dir} {input.assembly} {output}
    """
