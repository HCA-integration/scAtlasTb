#!/usr/bin/env bash

usage() {
  cat <<EOF
usage: $0 options

Install or update all environments in ./envs for the pipeline.

OPTIONS:
   -h     Show this message
   -c     Command to install conda packages, either 'mamba' or 'conda' (default: conda)
   -q     Quiet installation
   -n     Dry run (do not actually install/remove environments)
   -r     Remove all environments defined in ./envs and exit
EOF
}

CONDA_CMD="conda"
QUIET=""
ENVS_DIR=$(dirname "$0")
DRYRUN=""
REMOVE=0

while getopts "hc:nqr" OPTION; do
    case $OPTION in
        h) usage; exit 0;;
        c) CONDA_CMD=$OPTARG;;
        n) DRYRUN="-n";;
        q) QUIET="-q";;
        r) REMOVE=1;;
        ?) usage; exit 1;;
    esac
done

if [[ $REMOVE -eq 1 ]]; then
    echo "Removing environments defined in $ENVS_DIR:"
    for file in "$ENVS_DIR"/*.yaml; do
        if [[ -f "$file" ]]; then
            bash "$ENVS_DIR/remove_environment.sh" -f "$file" $DRYRUN -c "$CONDA_CMD"
        fi
    done
    $CONDA_CMD env list
    echo "Done removing all environments."
    exit 0
fi

if [[ $DRYRUN != "" ]]; then
    echo "Dry run: not installing environments..."
fi

for file in "$ENVS_DIR"/*.yaml; do 
    if [[ -f "$file" ]]; then
        bash "$ENVS_DIR/install_environment.sh" -f "$file" $DRYRUN -c "$CONDA_CMD" $QUIET
    fi 
done

# List all environments
$CONDA_CMD env list

echo "Done."
