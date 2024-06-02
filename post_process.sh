#!/bin/bash

log_date=$1
# assert that the log_date is not empty
if [ -z "$log_date" ]; then
    echo "log_date is empty"
    exit 1
fi

# generate results.csv and bestSol.csv
find /scratch/htc/spu/feasibilityjump/logs/"$log_date" -type f -name "*.log" -print0 | xargs -0 -n 1 awk -f fj_sol.awk

# post process of results.csv and bestSol.csv
# remove `\.opb.*\.log` from the first column of results.csv and bestSol.csv by awk
perl -pe 'BEGIN{$,=","} $_ = reverse $_; s/gol\..*?bpo\.//; $_ = reverse $_' results.csv > results.csv.tmp
mv results.csv.tmp results.csv
perl -pe 'BEGIN{$,=","} $_ = reverse $_; s/gol\..*?bpo\./bpo./; $_ = reverse $_' bestSol.csv > bestSol.csv.tmp
mv bestSol.csv.tmp bestSol.csv

# Read BIG-LIN.test and extract basenames without extension
while read -r line; do
    basename "$line" .opb
done < /scratch/htc/spu/github/roundingsat/BIG-LIN.test > big_lin_basenames.txt
# Remove duplicates
sort -u big_lin_basenames.txt -o big_lin_basenames.txt

# Read results.csv and extract basenames without extension
# awk -F',' '{gsub(/_.*$/, "", $1); print $1}' results.csv > results_basenames.txt
awk -F',' '{ print $1}' results.csv > results_basenames.txt
sort -u results_basenames.txt -o results_basenames.txt

# Compare the two lists and store the basenames that are in BIG-LIN.test but not in results.csv
not_run_basenames=$(grep -Fvxf results_basenames.txt big_lin_basenames.txt)

# Print the basenames
echo "$not_run_basenames"

# Count the number of basenames
not_run=$(echo "$not_run_basenames" | wc -l)
echo "Number of instances not run: $not_run"

# Write the basenames to a file
echo "$not_run_basenames" > not_run_basenames.txt

# Read not_run_basenames.txt
while read -r basename; do
    # Find the full path in BIG-LIN.test that ends with the basename
    grep "/${basename}.opb" /scratch/htc/spu/github/roundingsat/BIG-LIN.test
done < "$not_run_basenames" > BIG-LIN-NONSOLVED.test