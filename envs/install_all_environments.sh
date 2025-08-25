#!/usr/bin/env bash

usage() {
  cat <<EOF
usage: $0 options

Install or update all environments in ./envs for the pipeline.

OPTIONS:
   -h     Show this message
   -c     Command to install conda packages, either 'mamba' or 'conda' (default: conda)
   -q     Quiet installation
EOF
}

CONDA_CMD="conda"
QUIET=""
ENVS_DIR=$(dirname $0)
DRYRUN=""

while getopts "hc:nq" OPTION; do
    case $OPTION in
        h) usage; exit 1;;
        c) CONDA_CMD=$OPTARG;;
        n) DRYRUN="-n";;
        q) QUIET="-q";;
        ?) usage; exit;;
    esac
done

if [[ DRYRUN == "" ]]; then
    echo "Dry run: not installing environments..."
fi

for file in $ENVS_DIR/*; do 
    if [ -f "$file" ] && [[ $file == *yaml ]]; then
        bash $ENVS_DIR/install_environment.sh -f $file $DRYRUN -c $CONDA_CMD $QUIET
    fi 
done

# List all environments
$CONDA_CMD env list

echo "Done."
