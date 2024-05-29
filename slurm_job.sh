#!/bin/bash
#SBATCH --job-name=pbo_heur_test
#SBATCH --partition=opt_int
#SBATCH --time=0-00:31:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/htc/spu/feasibilityjump/build-Release/logs/%x_%A_%a.log
#SBATCH --error=/scratch/htc/spu/feasibilityjump/build-Release/errs/%x_%A_%a.log
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=56

set -x  # Enable printing of commands before execution

bin_dir=/scratch/htc/spu/feasibilityjump/build-Release
# Read each line from the file
while IFS= read -r line
do
  # Submit a job for each line
  srun --ntasks=1 --cpus-per-task=1 --job-name="$line" --mem=8000 -p opt_int -A optimi_integer --time=00-00:31:00 --cpu-freq=highm1 --output=${bin_dir}/logs/%x_%A_%a.log --error=${bin_dir}/errs/%x_%A_%a.log ${bin_dir}/xpress_fj -v -h -t 1800 "$line" &
done < "/scratch/htc/bzfmexig/pb_comp/OPT-LIN-BIG.test"