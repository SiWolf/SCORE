# Script for deleting existing SCORE results
# This is for testing various experimental combinations and parameters
# FastQC and Flexbar files should not need to be deleted
# To empty all results folders: ./empty_results.sh full
# To empty all results folders but save previous results: ./empty_results.sh full Test_Folder

path_parent_folder=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$path_parent_folder"

MODE=$1

cd simulation_data/
rm -r *.gz
rm -r *.idx
rm -r kallisto/sample*
rm -r *.rda
rm -r *.txt

cd ../../../deg/

if [ -z "$2" ]
then
	rm -r consensus*
	rm -r deg*
	rm -r diffexpr*
	rm -r filtered*
	rm -r pathway*
else	
	ANALYSIS=$2
	mkdir $ANALYSIS
	mv consensus* $ANALYSIS
	mv deg* $ANALYSIS
	mv diffexpr* $ANALYSIS
	mv filtered* $ANALYSIS
	mv pathway* $ANALYSIS
fi

if [ "$MODE" = "full" ]
then
	cd ../fastqc/
	rm -r *
fi

cd ../mapped/bowtie2/featureCounts/
rm -r *
cd ../../kallisto/
rm -r *
cd ../../trimmed/

if [ "$MODE" = "full" ]
then
	rm -r *
fi

cd ../
rm log.txt