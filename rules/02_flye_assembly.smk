
# Use flye for de-novo long reads assembly
rule flye_assembly:
  input: join(TMP, 'reads', '{strain}_long_reads_polished.fa'),
  output: assembly = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  params:
    flye_dir = directory(join(TMP, '{strain}', 'flye'))
  threads: 12
  singularity: "docker://cmdoret/flye:2.3.6"
  shell:
    """
    flye --nano-corrected {input} \
         --threads {threads} \
         --iterations 3 \
         -o {params.flye_dir} \
         -g 45m
    mv {params.flye_dir}/scaffolds.fa {output}
    """


# Use mummer to compute all-vs-all contig similarity
rule dnadiff_pairwise:
  input: join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa')
  output: join(TMP, '{strain}', 'contig_similarity.tsv')
  singularity: "docker://quay.io/biocontainers/mummer"
  script: "./scripts/pairwise_dnadiff {input} {output}"


# filter out small contigs of unmerged haplotypes
rule filter_haplotypes:
  input:
    genome = join(OUT, 'assemblies', '01_Ac_{strain}_flye.fa'),
    similarity = join(TMP, '{strain}', 'contig_similarity.tsv')
  output:
    keep = join(OUT, 'assemblies', '02_Ac_{strain}_haplofilter.fa'), 
    drop = join(TMP, '{strain}_drop_haplotypes.fa')
  params:
    similarity = 90
  run:
    sim = pd.read_csv(input['similarity'], sep='\t', header=None)
    sim.columns = ['contig1', 'contig2', 'aligned1', 'aligned2']
    sim['highest_sim'] = sim.apply(lambda x: max(x.aligned1, x.aligned2, axis=0))
    good_pairs = sim.loc[sim.highest_sim > params['similarity'],: ]
    fasta_keep = open(output['keep'])
    fasta_drop = open(output['drop'])
    with SeqIO.parse(input[0], 'fasta') as genome:
        for rec in genome:
            if rec.id in good_contigs:
                SeqIO.write(rec, fasta_keep, 'fasta')
            else:
                SeqIO.write(rec, fasta_drop, 'fasta')
        {params.similarity} \

