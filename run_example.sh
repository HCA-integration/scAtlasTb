#!/usr/bin/env bash
set -e -x

snakemake \
  --profile .profiles/local \
  --configfile configs/quickstart.yaml \
  --snakefile workflow/Snakefile \
    $@
