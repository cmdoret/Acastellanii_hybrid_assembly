# This is an optional workflow used to separate mitochondrial from nuclear contigs prior to scaffolding.
# This is done to prevent the mitochondrial genome from being incorporated into another contig.

# Align contigs to the mitochondrial genome from previous assembly 
# to find which are mitochondrial
rule split_mito_tigs:
    input:
        genome = join(OUT, 'assemblies', '02_Ac_{strain}_hypo.fa'),
        mito= config['mitochondrion'] 
    output:
        nucl= join(OUT, 'assemblies', '02_Ac_{strain}_hypo_nucl.fa'),
        mito = join(OUT, 'assemblies', '02_Ac_{strain}_hypo_mito.fa')

    conda: "../envs/align.yaml"
    shell:
        """
        python scripts/split_mito.py \
            -g {input.genome} \
            -m {input.mito} \
            -N {output.nucl} \
            -M {output.mito}
        """
