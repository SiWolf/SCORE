# ------------------------------
# Title: empty_results.sh
# Author: Silver A. Wolf
# Last Modified: Mo, 09.12.2019
# Version: 0.0.1
# ------------------------------

# Script for deleting existing SCORE results
# This is used for testing various experimental combinations and parameters
# To empty all primary results: ./empty_results.sh default
# To empty all results: ./empty_results.sh full
# To empty all primary results but save previous results: ./empty_results.sh default Test_Folder

PARENT_FOLDER_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
MODE=$1

# Simulation data is completely removed
cd "$PARENT_FOLDER_PATH"
cd simulation_data/
rm -r *.gz
rm -r *.idx
rm -r kallisto/sample*
rm -r *.rda
rm -r *.txt
cd ../../../deg/

# DEG results are either deleted or moved to a specified directory
if [ -z "$2" ]
then
	rm -r consensus*
	rm -r dag*
	rm -r deg*
	rm -r diffexpr*
	rm -r filtered*
	rm -r genes*
	rm -r hmmer*
	rm -r log*
	rm -r pathway*
	rm -r proteins*
	rm -r transcript*
	rm -r summary*
else	
	ANALYSIS=$2
	mkdir $ANALYSIS
	mv consensus* $ANALYSIS
	mv dag* $ANALYSIS
	mv deg* $ANALYSIS
	mv diffexpr* $ANALYSIS
	mv filtered* $ANALYSIS
	mv genes* $ANALYSIS
	mv hmmer* $ANALYSIS
	mv log* $ANALYSIS
	mv pathway* $ANALYSIS
	mv proteins* $ANALYSIS
	mv summary* $ANALYSIS
	mv transcript* $ANALYSIS
fi

rm -r benchmarking*
cd ../
rm log.txt
rm references/*.tmp

# FastQC, flexbar and MultiQC files are deleted on default (default mode)
if [ "$MODE" = "default" ]
then
	cd fastqc/
	rm -r *
	cd ../trimmed/
	rm -r *
	cd ../
fi

# In addition, Bowtie2, featureCounts and kallisto files are deleted on full reset (full mode)
if [ "$MODE" = "full" ]
then
	cd fastqc/
	rm -r *
	cd ../trimmed/
	rm -r *
	cd ../mapped/bowtie2/
	rm -r *.sam
	cd featureCounts/
	rm -r *
	cd ../../kallisto/
	rm -r *
	cd ../../
fi