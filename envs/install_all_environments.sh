#!/usr/bin/env bash

usage() {
  cat <<EOF
usage: $0 options

Install or update all environments in ./envs for the pipeline.

OPTIONS:
   -h     Show this message
   -m     Command to install conda packages, either 'mamba' or 'conda' (default: mamba)
   -q     Quiet installation
EOF
}

MAMBA_CMD="mamba"
QUIET=""
ENVS_DIR=$(dirname $0)
DRYRUN=""

while getopts "hm:nq" OPTION; do
    case $OPTION in
        h) usage; exit 1;;
        m) MAMBA_CMD=$OPTARG;;
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
        bash $ENVS_DIR/install_environment.sh -f $file $DRYRUN -m $MAMBA_CMD $QUIET
    fi 
done

# List all environments
$MAMBA_CMD env list

echo "Done."
