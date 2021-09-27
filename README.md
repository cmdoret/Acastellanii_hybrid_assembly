### Acanthamoeba castellanii Hybrid assembly pipeline

This repository contains an automatic and pipeline for the genome assembly of Acanthamoeba castellanii. The different steps of the pipeline are illustrated below.
The analysis is implemented by combining the snakemake workflow system with [singularity containers](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#singularity) and conda environments.

#### Dependencies

 * python >=3.7
 * singularity >= 3.2
 * conda >= 4.8

Python packages:

 * snakemake >= 5.5

#### Usage

The pipeline can be started using:

```bash
snakemake --use-conda -j4
```
or

```bash
snakemake --use-singularity -j4
```
to use singularity instead of conda environments.

#### Pipeline steps

The pipeline uses 3 types of input data:
 * shotgun Illumina reads
 * Hi-C Illumina reads
 * Oxford Nanopore long reads

All those reads are automatically downloaded from SRA when running the pipeline.

The initial assembly is performed with long reads only using Flye. The short reads are then used to polish this assembly using HyPo. The Hi-C scaffolding is done using instagraal, followed by instagraal-polish to fix errors introduced by instagraal.

>TODO: Add quast report at the end of the pipeline

Each rule requiring a third party software pulls a standalone container hosted on dockerhub or quay to work in an isolated environment with a fixed version of the software.

![image](doc/assembly.svg)

Unfortunately, instagraal requires access to a GPU with CUDA drivers. It is currently not possible to make it compatible with singularity. This means instagraal has to be installed on the host machine for the scaffolding to work.

### References
Tools used in this pipeline:
 * filtlong v0.2.0: 
 * [flye](https://github.com/fenderglass/Flye/) v2.3.6: doi:10.1073/pnas.1604560113
 * [seqtk](https://github.com/lh3/seqtk) v1.3
 * [bowtie2](https://github.com/BenLangmead/bowtie2) v2.3.4.1: doi:10.1038/nmeth.1923
 * [minimap2](https://github.com/lh3/minimap2) v2.17: doi:10.1093/bioinformatics/bty191
 * [hypo](https://github.com/kensung-lab/hypo) v0.1.0: doi:10.1101/2019.12.19.882506
 * [instagraal](https://github.com/koszullab/instaGRAAL) v0.1.6: doi:10.1038/ncomms6695
 * [hicstuff](https://github.com/koszullab/hicstuff) v2.2.2: doi:10.5281/zenodo.2620608
