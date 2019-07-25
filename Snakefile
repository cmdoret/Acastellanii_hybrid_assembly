# Docker-based pipeline for Acanthamoeba castellanii genome assembly orchestrated using snakemake
# cmdoret, 20190718

from os.path import join
from snakemake.utils import validate
import pandas as pd


## CONFIGURATION FILES ##


configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("strain", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["strain", "unit"], drop=False)
# Enforces str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

TMP = config['tmp_dir']
OUT = config['out_dir']

## WILDCARD CONSTRAINTS
wildcard_constraints:
  sample="|".join(samples.index)

## Helper functions

def get_fastqs(wildcards, libtype="shotgun"):
    """Get fastq files of a particular library type from unit sheet"""
    fqs = units.loc[(units.strain == wildcards.strain) & (units.libtype == libtype), ["fq1", "fq2"]].dropna()
    if len(fqs.values) == 2:
        return {"r1": fqs.fq1.values, "r2": fqs.fq2.values}
    return {"r1": fqs.fq1.values}

## PIPELINE
include: "rules/01_reads_processing.smk"
include: "rules/02_flye_assembly.smk"
include: "rules/03_pilon_racon_polishing.smk"
include: "rules/04_instagraal_scaffolding.smk"
include: "rules/05_pilon_polishing2.smk"
# include: "rules/06_LINKS_scaffolding.smk"
# include: "rules/07_quast_report.smk"

rule all:
    input: expand(join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa'), strain=samples['strain'])

