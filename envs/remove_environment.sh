#!/usr/bin/env bash
set -e

usage() {
  cat <<EOF
usage: $0 options

Remove an environment defined by a YAML file or by environment name

OPTIONS:
   -h     Show this message
   -f     Environment YAML file OR environment name
   -c     Command to manage conda environments, either 'mamba' or 'conda' (default: conda)
   -n     Dry run (do not actually remove environment)
EOF
}

CONDA_CMD="conda"
TARGET=""
EXECUTE=true

while getopts "hf:c:n" OPTION; do
    case $OPTION in
        h) usage; exit 0;;
        f) TARGET=$OPTARG;;
        c) CONDA_CMD=$OPTARG;;
        n) EXECUTE=false;;
        ?) usage; exit 1;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    esac
done

if [[ -z "$TARGET" ]]; then
    echo "Error: Must specify environment YAML file or environment name with -f" >&2
    usage
    exit 1
fi

# If TARGET is a file, extract environment name from YAML
if [[ -f "$TARGET" ]]; then
    ENV_NAME=$(grep -E '^name:' "$TARGET" | awk '{print $2}')
else
    ENV_NAME="$TARGET"
fi

if ! $CONDA_CMD env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
    echo "Environment $ENV_NAME does not exist, skipping"
    exit 0
fi

if ! $EXECUTE; then
    echo "[DRYRUN] Would remove environment: $ENV_NAME"
else
    echo "Removing environment: $ENV_NAME"
    $CONDA_CMD env remove -n "$ENV_NAME" -y
fi
