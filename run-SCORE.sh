# Run script for the Snakemake file
# Contains basic settings for running the Snakefile
# The first given argument is the number of cores Snakemake will be able to use
# Default number of cores used is 4
# How to run: ./SCORE.sh <#cores>

if [ -z "$1" ]
then
	time snakemake -s "SCORE.snk" -j 4 --use-conda
else
	time snakemake -s "SCORE.snk" -j $1 --use-conda
fi