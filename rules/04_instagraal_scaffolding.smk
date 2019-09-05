
# Generate Hi-C matrix from raw reads
rule hicstuff_hic_processing:
    input:
      r1 = GS.remote(join(TMP, "reads", "{strain}_Hi-C.end1.fq.gz")),
      r2 = GS.remote(join(TMP, "reads", "{strain}_Hi-C.end2.fq.gz")),
      assembly = GS.remote(join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa'))
    output: GS.remote(directory(join(TMP, "hicstuff", "{strain}")))
    threads: CPUS
    resources: mem=32000
    params:
      enzyme = "DpnII",
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
    assembly = GS.remote(join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')),
    hicstuff_dir = GS.remote(join(TMP, "hicstuff", "{strain}"))
  output: GS.remote(join(OUT, 'assemblies', '05_Ac_{strain}_instagraal.fa'))
  params:
    instagraal_input_dir = GS.remote(join(TMP, 'instagraal', '{strain}'))
  shell:
    """
    mkdir {params.instagraal_input_dir}
    cp "{input.hicstuff_dir}/abs_fragments_contacts_weighted.txt" \
       "{input.hicstuff_dir}/fragments_list.txt" \
       "{input.hicstuff_dir}/info_contigs.txt" \
       "{params.instagraal_input_dir}"
    instagraal {params.instagraal_input_dir} {input.assembly} {output}
    """
