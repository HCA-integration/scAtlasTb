#!/usr/bin/env bash
set -e -x

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/metrics/hlcav1.yaml \
  --snakefile workflow/Snakefile \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
    $@
