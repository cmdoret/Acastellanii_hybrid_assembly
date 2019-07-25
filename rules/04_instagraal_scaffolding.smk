
# Generate Hi-C matrix from raw reads
rule hicstuff_hic_processing:
    input:
      unpack(lambda w: get_fastqs(w, libtype='Hi-C')),
      assembly = join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa')
    output: directory(join(TMP, "hicstuff", "{strain}"))
    threads: 12
    params:
      enzyme = "DpnII"
    shell:
      """
      hicstuff pipeline -t {threads} \
                        -e {params.enzyme} \
                        -g {input.assembly} \
                        {input.hic_ends} \
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
