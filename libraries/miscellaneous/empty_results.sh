# -------------------------------
# Title: empty_results.sh
# Author: Silver A. Wolf
# Last Modified: Wed, 18.12.2019
# Version: 0.0.3
# -------------------------------

# Script for deleting existing SCORE results
# Useful for rapidly testing various experimental combinations and parameters
# Empty all/none/primary results: ./empty_results.sh <full/fair/default> 
# Empty all/none/primary results and save previous run: ./empty_results.sh <full/fair/default> <Analysis_Folder>

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
	ANALYSIS=$2-`date '+%Y-%m-%d-%H-%M'`-`uuidgen -t | head -c 5`
	
	# Iff an older analysis folder is found, remove it
	if [ -d "$ANALYSIS" ]
	then
		rm -r $ANALYSIS
	fi
	
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
rm tmp.db
rm log.txt
rm references/*.fai
rm references/*.tmp

# Fair mode
# Does not delete any of the results

# Default mode
# FastQC, flexbar and MultiQC files are deleted on default
if [ "$MODE" = "default" ]
then
	cd fastqc/
	rm -r *
	cd ../trimmed/
	rm -r *
	cd ../
fi

# Full mode
# In addition, Bowtie2, featureCounts and kallisto files are deleted on full reset
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