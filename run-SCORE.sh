# -----------------------------
# Title: run-SCORE.sh
# Author: Silver A. Wolf
# Last Modified: Fr, 06.12.2019
# Version: 0.0.1
# -----------------------------

# Run script for the Snakemake file
# Contains basic settings for running the Snakefile
# The first given argument is the number of cores Snakemake will be able to use
# Default number of cores used is 4
# How to run: ./SCORE.sh <#cores>
# Will automatically empty previous result-files in order to ensure no collisions

./libraries/miscellaneous/empty_results.sh full

snakemake -s "SCORE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg

if [ -z "$1" ]
then
	{ time snakemake -s "SCORE.snk" -j 4 --use-conda ; } 2>&1 | tee deg/log.txt
else
	{ time snakemake -s "SCORE.snk" -j $1 --use-conda ; } 2>&1 | tee deg/log.txt
fi