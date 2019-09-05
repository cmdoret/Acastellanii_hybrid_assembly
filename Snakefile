# Docker-based pipeline for Acanthamoeba castellanii genome assembly orchestrated using snakemake
# cmdoret, 20190718

from os.path import join
from snakemake.utils import validate
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
import pandas as pd
import numpy as np

## CONFIGURATION FILES ##

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep='\t').set_index("strain", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep='\t', dtype=str).set_index(["strain", "unit"], drop=False)
# Enforces str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

bucket = config['remote']['bucket']
GS = GSRemoteProvider()
CPUS = config['n_cpus']
TMP = join(bucket, config['tmp_dir'])
OUT = join(bucket, config['out_dir'])

## WILDCARD CONSTRAINTS
wildcard_constraints:
  strain="|".join(samples.index),
  libtype="|".join(np.unique(units.libtype))

## Helper functions

def access_remote(local_path):
    """
    Given the local path to a file, return the full path
    for remote access by accessing snakemake config. When using
    SFTP (ssh access to remote files), the bucket parameter is used
    as a prefix to the path.

    Parameters
    ----------
    local_path : str
        The local path within filesystem that need to be wrapped.

    Returns
    -------
        The path with remote information (provider, host, bucket name) prepended.

    """
    rc = config['remote']
    provider = rc['provider']
    bucket = rc['bucket']
    host = rc['host']
    username = rc['username']
    ssh_key = rc['ssh_key']
    password = rc['password']

    if not isinstance(local_path, str):
      print("Remote information cannot be added to expanded path. Expand after access_remote.")
      sys.exit(1)

    if provider == "GS":
      GS = GSRemoteProvider()
      remote_path = GS.remote(join(bucket, local_path))
    elif provider == 'SFTP':
      if password == "":
        SFTP = SFTPRemoteProvider(username=username, private_key=ssh_key)
      else:
        SFTP = SFTPRemoteProvider(username=username, password=password)
      remote_path = SFTP.remote(host + ":22" + join(bucket, local_path))
    elif provider in ('', "local"):
      remote_path = local_path
    else:
      print('Pipeline not configured for remote provider: {}'.format(provider))
      print('Please edit the config.yaml file to select a valid provider.')
    return remote_path


def get_fastqs(wildcards):
    """
    Get fastq files (units) of a particular library type of one sample 
    from the unit sheet
    """
    fqs = units.loc[(units.strain == wildcards.strain) & (units.libtype == wildcards.libtype), ["fq1", "fq2"]]
    fq1 = fqs.fq1.dropna()
    fq2 = fqs.fq2.dropna()
    if len(fq2) == len(fq1):
        return {
          "r1": list(map(access_remote, fqs.fq1)),
          "r2": list(map(access_remote, fqs.fq2))
        }
    return {"r1": list(map(access_remote, fqs.fq1))}

def ok_google():
  """
  Used to work around the incompatibility between google cloud and functions as input.
  Builds a dictionary of paths for librairies and libtypes and prepends google storage 
  bucket name to each path.
  """
  # path_dict = {libtype: {library: {fq1: [path1, ...], fq2: [path1, ...]}}}
  path_dict = {}
  for t in np.unique(units.libtype):
    path_dict[t] = {}
    for s in samples.strain:
      path_dict[t][s] = {}
      for e in ["fq1", "fq2"]:
        fqs = units.loc[(units.strain == s) & (units.libtype == t), e].dropna().values.tolist()
        path_dict[t][s][e] = list(map(lambda x: join(bucket, x), fqs))
  return path_dict

units_dict = ok_google()

## PIPELINE
include: "rules/01_reads_processing.smk"
include: "rules/02_flye_assembly.smk"
include: "rules/03_pilon_racon_polishing.smk"
include: "rules/04_instagraal_scaffolding.smk"
include: "rules/05_pilon_polishing2.smk"
# include: "rules/06_LINKS_scaffolding.smk"
# include: "rules/07_quast_report.smk"

rule all:
    input: GS.remote(expand(join(OUT, 'assemblies', '06_Ac_{strain}_pilon2.fa'), strain=samples['strain']))
