#!/usr/bin/env bash
set -e -x

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/metrics/quickstart.yaml \
  --snakefile workflow/Snakefile \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
    $@
