#!/bin/bash
#SBATCH --job-name=glyverity_census
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
# Census of one run as a batch job on Sapelo2 (GACRC).
# Usage (from the repository folder):
#   sbatch scripts/slurm/census.sbatch <run.d or run.mzML> <output folder> [census_run options]
# Results go to the output folder; the printed report goes to <job name>.<job id>.out and
# stage times to the .err file, both in the folder where sbatch was called.
set -euo pipefail
source "$HOME/venvs/lcms/bin/activate"
cd "$SLURM_SUBMIT_DIR"
echo "job $SLURM_JOB_ID on $(hostname), $(date)"
echo "glyverity commit $(git rev-parse --short HEAD)"
python -c "import sys, numpy; print('python', sys.executable, 'numpy', numpy.__version__, numpy.__file__)"
mkdir -p "$2"
python -u scripts/census_run.py "$@"
echo "finished $(date)"
