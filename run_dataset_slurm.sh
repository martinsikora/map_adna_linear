#!/usr/bin/env bash

set -euo pipefail

# Usage information
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
${COLOR_BOLD}${COLOR_BLUE}Usage:${COLOR_RESET} $(basename "$0") ${COLOR_GREEN}DATASET${COLOR_RESET} [OPTIONS]

Submit Snakemake workflow to SLURM with automatic resource management.

${COLOR_BOLD}${COLOR_BLUE}Arguments:${COLOR_RESET}
  ${COLOR_GREEN}DATASET${COLOR_RESET}       Path to dataset directory containing config/config.yml

${COLOR_BOLD}${COLOR_BLUE}Options:${COLOR_RESET}
  ${COLOR_YELLOW}-j, --jobs N${COLOR_RESET}          Maximum number of jobs to run in parallel (default: 100)
  ${COLOR_YELLOW}-p, --partition NAME${COLOR_RESET}  SLURM partition to use (default: general)
  ${COLOR_YELLOW}-t, --time TIME${COLOR_RESET}       Max walltime per job (default: 4:00:00)
  ${COLOR_YELLOW}-m, --mem MB${COLOR_RESET}          Default memory per job in MB (default: 8000)
  ${COLOR_YELLOW}-c, --config FILE${COLOR_RESET}     Alternate config file (default: config/config.yml)
  ${COLOR_YELLOW}--dry-run${COLOR_RESET}             Show what would be done without submitting jobs
  ${COLOR_YELLOW}--unlock${COLOR_RESET}              Unlock working directory
  ${COLOR_YELLOW}--${COLOR_RESET}                   Pass remaining args directly to Snakemake
  ${COLOR_YELLOW}-h, --help${COLOR_RESET}            Show this help message

${COLOR_BOLD}${COLOR_BLUE}Example:${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example --jobs 50 --partition highmem${COLOR_RESET}

  # For testing without submitting:
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example --dry-run${COLOR_RESET}
  ${COLOR_EXAMPLE}bash $(basename "$0") dataset_example -- --rerun-incomplete --keep-going${COLOR_RESET}

EOF
    exit 0
}

# Default parameters
MAX_JOBS=100
PARTITION="general"
MAX_TIME="4:00:00"
DEFAULT_MEM=8000
SNAKEMAKE_ARGS=()
UNLOCK=""
CONFIG_FILE="config/config.yml"

# Parse command line arguments
if [[ $# -eq 0 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
fi

DATASET="$1"
shift

while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jobs)
            MAX_JOBS="$2"
            shift 2
            ;;
        -p|--partition)
            PARTITION="$2"
            shift 2
            ;;
        -t|--time)
            MAX_TIME="$2"
            shift 2
            ;;
        -m|--mem)
            DEFAULT_MEM="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --dry-run)
            SNAKEMAKE_ARGS+=(--dry-run)
            shift
            ;;
        --unlock)
            UNLOCK="--unlock"
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
            if [[ "$1" == -* ]]; then
                SNAKEMAKE_ARGS+=("$1")
                if [[ $# -ge 2 ]] && [[ "$2" != -* ]]; then
                    SNAKEMAKE_ARGS+=("$2")
                    shift 2
                else
                    shift
                fi
            else
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
            fi
            ;;
    esac
done

# Get the directory where this script is located
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Setup paths - workflow is relative to script location
WORKFLOW_DIR="$SCRIPT_DIR/workflow"
RUN_DIR="$(realpath $DATASET)"

# Change to the run directory
cd "$RUN_DIR"

echo "=========================================="
echo "Snakemake SLURM Submission"
echo "=========================================="
echo "Dataset:       $RUN_DIR"
echo "Max jobs:      $MAX_JOBS"
echo "Partition:     $PARTITION"
echo "Default time:  $MAX_TIME"
echo "Default mem:   ${DEFAULT_MEM}M"
echo "=========================================="

# Run unlock if requested
if [[ -n "$UNLOCK" ]]; then
    echo "Unlocking working directory..."
    snakemake \
      -s "$WORKFLOW_DIR/Snakefile" \
      --configfile "$CONFIG_FILE" \
      --unlock
    echo "Directory unlocked."
    exit 0
fi

# Convert time format (HH:MM:SS to minutes for runtime parameter)
RUNTIME_MIN=$(echo "$MAX_TIME" | awk -F: '{ print ($1 * 60) + $2 }')

# Run Snakemake with SLURM executor (Snakemake 9+)
snakemake \
  -s "$WORKFLOW_DIR/Snakefile" \
  --configfile "$CONFIG_FILE" \
  --executor slurm \
  --jobs ${MAX_JOBS} \
  --default-resources slurm_partition=${PARTITION} mem_mb=${DEFAULT_MEM} runtime=${RUNTIME_MIN} \
  --retries 2 \
  --latency-wait 60 \
  --printshellcmds \
  "${SNAKEMAKE_ARGS[@]}"

echo "=========================================="
echo "Workflow submission complete!"
echo "Monitor jobs with: squeue -u \$USER"
echo "Workflow logs: $RUN_DIR/logs/"
echo "=========================================="
