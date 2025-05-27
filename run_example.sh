#!/usr/bin/env bash
set -e -x

snakemake \
  --profile .profiles/local \
  --configfile \
    configs/example_config.yaml \
  --snakefile workflow/Snakefile \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
    $@
