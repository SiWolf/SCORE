# Run script for the Snakemake file
# Contains basic settings for running the Snakefile
# The first given argument is the number of cores Snakemake will be able to use
# Default number of cores used is 4
# How to run: ./SCORE.sh <#cores>
# Will automatically empty previous result-files in order to ensure no collisions

./libraries/miscellaneous/empty_results.sh full

if [ -z "$1" ]
then
	{ time snakemake -s "SCORE.snk" -j 4 --use-conda ; } 2>&1 | tee log.txt
else
	{ time snakemake -s "SCORE.snk" -j $1 --use-conda ; } 2>&1 | tee log.txt
fi
