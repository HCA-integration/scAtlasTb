#!/usr/bin/env bash
set -e

usage() {
  cat <<EOF
usage: $0 options

Install or update an environment

OPTIONS:
   -h     Show this message
   -f     Environment YAML file
   -m     Command to install conda packages, either 'mamba' or 'conda' (default: mamba)
   -q     Quiet installation
EOF
}

MAMBA_CMD="mamba"
QUIET="" # "-q"
ENVS_DIR=$(dirname $0)
EXECUTE=true

while getopts "hf:m:nq" OPTION; do
    case $OPTION in
        h) usage; exit 1;;
        f) FILE=$OPTARG;;
        m) MAMBA_CMD=$OPTARG;;
        n) EXECUTE=false;;
        q) QUIET="-q";;
        ?) usage; exit;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    esac
done

if [[ $FILE == "" ]]; then
    echo "Error: Must specify environment YAML file with -f" >&2
    usage

    exit 1
fi

ENV=$(basename $FILE | cut -d. -f1)
ALL_ENVS=$($MAMBA_CMD env list)

if [[ $ALL_ENVS == *"$ENV"* ]]; then
    operation="update"
else
    operation="create"
fi
echo "$operation $ENV from $FILE..."
if $EXECUTE; then
    $MAMBA_CMD env $operation $QUIET -f $FILE
fi