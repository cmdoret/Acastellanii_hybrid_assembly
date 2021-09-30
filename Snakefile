# Docker-based pipeline for Acanthamoeba castellanii genome assembly orchestrated using snakemake
# cmdoret, 20190718

from os.path import join
from snakemake.utils import validate
import pandas as pd
import numpy as np

## CONFIGURATION FILES ##

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t', comment='#').set_index("strain", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep='\t', dtype=str).set_index(["strain", "unit"], drop=False)
# Enforces str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

CPUS = config['n_cpus']
TMP = config['tmp_dir']
OUT = config['out_dir']

## WILDCARD CONSTRAINTS
wildcard_constraints:
  strain="|".join(samples.index),
  libtype="|".join(np.unique(units.libtype))

## Helper functions


conda: "envs/genome_assembly.yaml"
## PIPELINE
include: "rules/00_downloaders.smk"
include: "rules/01_reads_processing.smk"
include: "rules/02_flye_assembly.smk"
include: "rules/03_polishing.smk"
include: "rules/03a_mitochondria.smk"
include: "rules/04_instagraal_scaffolding.smk"
include: "rules/05_pilon_polishing2.smk"
# include: "rules/06_LINKS_scaffolding.smk"
# include: "rules/07_quast_report.smk"

rule all:
  input:
    expand(join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa'), strain=samples['strain'])
    #expand(join(OUT, 'assemblies', '04_Ac_{strain}_instagraal_polish.fa'), strain=samples['strain'])


