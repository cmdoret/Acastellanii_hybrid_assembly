# This is an optional workflow used to separate mitochondrial from nuclear contigs prior to scaffolding.
# This is done to prevent the mitochondrial genome from being incorporated into another contig.

# Align contigs to the mitochondrial genome from previous assembly 
# to find which are mitochondrial
rule split_mito_tigs:
    input:
        genome = join(OUT, 'assemblies', '04_Ac_{strain}_racon.fa'),
        mitochondrion = config['mitochondrion'] 
    output:
        nuclear = join(OUT, 'assemblies', '04a_Ac_{strain}_racon_nuclear.fa'),
        mito = join(OUT, 'assemblies', '04b_Ac_{strain}_racon_mito.fa')

    shell:
        """
        python scripts/split_mito.py \
            -g {input.genome} \
            -m {input.mitochondrion} \
            -N {output.nuclear} \
            -M {output.mito}
        """
