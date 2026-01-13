#!/usr/bin/env bash
set -euo pipefail

if [[ -t 1 ]]; then
    COLOR_BOLD="$(printf '\033[1m')"
    COLOR_EXAMPLE="$(printf '\033[90m')"
    COLOR_BLUE="$(printf '\033[34m')"
    COLOR_GREEN="$(printf '\033[32m')"
    COLOR_YELLOW="$(printf '\033[33m')"
    COLOR_RESET="$(printf '\033[0m')"
else
    COLOR_BOLD=""
    COLOR_EXAMPLE=""
    COLOR_BLUE=""
    COLOR_GREEN=""
    COLOR_YELLOW=""
    COLOR_RESET=""
fi

usage() {
    cat << EOF
${COLOR_BOLD}${COLOR_BLUE}Usage:${COLOR_RESET} $(basename "$0") ${COLOR_GREEN}DATASET${COLOR_RESET} [OPTIONS] [-- ${COLOR_YELLOW}SNAKEMAKE_ARGS${COLOR_RESET}...]

Run the Snakemake workflow locally.

${COLOR_BOLD}${COLOR_BLUE}Arguments:${COLOR_RESET}
  ${COLOR_GREEN}DATASET${COLOR_RESET}       Path to dataset directory containing config/config.yml

${COLOR_BOLD}${COLOR_BLUE}Options:${COLOR_RESET}
  ${COLOR_YELLOW}--dry-run${COLOR_RESET}     Run Snakemake with --dry-run
  ${COLOR_YELLOW}-c, --config FILE${COLOR_RESET}  Alternate config file (default: config/config.yml)
  ${COLOR_YELLOW}-h, --help${COLOR_RESET}    Show this help message

${COLOR_BOLD}${COLOR_BLUE}Examples:${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example -- --cores 8 --rerun-incomplete${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example --dry-run${COLOR_RESET}
EOF
    exit 0
}

if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
fi

DATASET="$1"
shift  # Remove the first argument (DATASET) from the argument list

SNAKEMAKE_ARGS=()
CONFIG_FILE="config/config.yml"

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --dry-run)
            SNAKEMAKE_ARGS+=(--dry-run)
            shift
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                SNAKEMAKE_ARGS+=("$1")
                shift
            done
            ;;
        *)
            SNAKEMAKE_ARGS+=("$1")
            shift
            ;;
    esac
done

WORKFLOW_DIR="$(realpath workflow)"
RUN_DIR="$(realpath "$DATASET")"

cd "$RUN_DIR"

snakemake \
  -s "$WORKFLOW_DIR/Snakefile" \
  --configfile "$CONFIG_FILE" \
  "${SNAKEMAKE_ARGS[@]}"
