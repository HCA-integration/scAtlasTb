#!/usr/bin/env bash
set -e -x

WORKDIR=$(dirname $(dirname $0))
cd $WORKDIR

#--snakefile $WORKDIR/Snakefile
snakemake --rerun-triggers mtime params input code --rerun-incomplete --configfile test/config.yaml --use-conda $@

conda run -n scanpy python test/run_assertions.py