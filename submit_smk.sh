#!/usr/bin/env bash
#SBATCH --no-requeue
#SBATCH --mem 6G
#SBATCH -p genoa64
#SBATCH --qos pipelines
#SBATCH --output=logs/slurm/smk.%j.out
#SBATCH --output=logs/slurm/smk.%j.err

# Configure bash
set -e          # exit immediately on error
set -u          # exit immidiately if using undefined variables
set -o pipefail # ensure bash pipelines return non-zero status if any of their command fails

# Setup trap function to be run when canceling the pipeline job. It will propagate the SIGTERM signal
# to snakemake so that all jobs launche by the pipeline will be cancelled too.
_term() {
        echo "Caught SIGTERM signal!"
        kill -s SIGTERM $pid
        wait $pid
}

trap _term TERM

# environment
eval "$(conda shell.bash hook)"
conda activate snakemake # just an env with snakemake (>8)

# Run the pipeline. The command uses the arguments passed to this script, e.g:
#
# $ sbatch submit_smk.sh --profile profiles/slurm --snakefile workflow/rules/explore.smk --forceall
#
# will use "--profile profiles/slurm workflow/rules/explore.smk --forceall" as arguments
snakemake --executor slurm "$@" & pid=$!

# Wait for the pipeline to finish
echo "Waiting for ${pid}"
wait $pid

# Return 0 exit-status if everything went well
exit 0

