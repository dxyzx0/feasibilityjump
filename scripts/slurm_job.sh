#!/bin/bash

set -x  # Enable printing of commands before execution

root_dir=/scratch/htc/spu/feasibilityjump
# get the current date and time
now=$(date +"%Y%m%d_%H%M%S")
# Read each line from the file
while IFS= read -r line
do
  # Submit a job for each line
  srun --exclusive --ntasks=1 --cpus-per-task=8 --job-name="$(basename "$line")" --mem=94000 -p opt_int -A optimi_integer --time=00-00:31:00 --cpu-freq=highm1 --output="${root_dir}/logs/${now}/%x_%A_%a.log" --error="${root_dir}/logs/${now}/%x_%A_%a.err" ${root_dir}/build-Release/pbo_fj -v -h -t 1800 /scratch/htc/spu/github/roundingsat/PB16-BIG/"${line}" &
done < "/scratch/htc/spu/github/roundingsat/BIG-LIN.test"
