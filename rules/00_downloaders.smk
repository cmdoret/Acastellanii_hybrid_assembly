# Rules for downloading data from the web

rule sra_dl_fq:
  message: "Getting {params.acc} into {output}"
  output: join('fq', '{strain}','{libtype}', '{libname}')
  params:
    acc = lambda w: units.sra[units.fq1.str.contains(w.libname)].values[0]
  conda: '../envs/sra.yaml'
  singularity: 'quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0'
  threads: 4
  shell:
    """
    fq={output}
    trim=${{fq%_[12].fastq.gz}}
    echo "SRA download to ${{trim}}_1.fastq and ${{trim}}_2.fastq"
    fasterq-dump -e {threads} "{params.acc}" -o $trim
    gzip ${{trim}}*fastq
    """