# -----------------------------
# Title: run-SCORE.sh
# Author: Silver A. Wolf
# Last Modified: Wed, 19.02.2020
# Version: 0.0.4
# -----------------------------

# Run script for SCORE
# Contains basic settings for running the Snakefile
# First given argument is the number of cores Snakemake will be able to use
# Second argument is SE or PE flag
# Will automatically empty previous result-files in order to ensure no collisions
# Default number of cores: 4
# How to run: ./run-SCORE.sh <#cores> <PE/SE>

./libraries/miscellaneous/empty_results.sh full

if [ -z "$1" ]
then
	snakemake -s "SCORE_PE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
	{ time snakemake -s "SCORE_PE.snk" -j 4 --use-conda ; } 2>&1 | tee deg/log.txt
else
	snakemake -s "SCORE_"$2".snk" --forceall --dag | dot -Tsvg > deg/dag.svg
	{ time snakemake -s "SCORE_"$2".snk" -j $1 --use-conda ; } 2>&1 | tee deg/log.txt
fi