#!/usr/bin/env bash
# Installs only the environments required to run configs/quickstart.yaml
# via run_example.sh (preprocessing_all integration_all metrics_all, use_gpu: false).
# Much faster than install_all_environments.sh, which installs all 20+ environments.
set -e

usage() {
  cat <<EOF
usage: $0 options

Install the conda environments required for the quickstart workflow used in run_example.sh.

OPTIONS:
   -h     Show this message
   -c     Command to install conda packages, either 'mamba' or 'conda' (default: conda)
   -q     Quiet installation
EOF
}

CONDA_CMD="conda"
QUIET=""
ENVS_DIR=$(dirname "$0")

while getopts "hc:q" OPTION; do
    case $OPTION in
        h) usage; exit 0;;
        c) CONDA_CMD=$OPTARG;;
        q) QUIET="-q";;
        ?) usage; exit 1;;
    esac
done

REQUIRED_ENVS=(
    snakemake
    scanpy
    bbknn
    scanorama
    scvi-tools
    scib
    pegasus
    plots
    funkyheatmap
)

for env in "${REQUIRED_ENVS[@]}"; do
    bash "$ENVS_DIR/install_environment.sh" -f "$ENVS_DIR/${env}.yaml" -c "$CONDA_CMD" $QUIET
done

$CONDA_CMD env list
echo "Done. Quickstart environments installed."
