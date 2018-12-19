# Script for deleting existing SCORE results
# This is for testing various experimental combinations and parameters
# FastQC and Flexbar files should not need to be deleted

cd simulation_data/
rm -r *.gz
rm -r *.idx
rm -r kallisto/sample*
rm -r *.rda
rm -r *.txt
cd ../../../deg/
rm -r consensus*
rm -r deg*
rm -r diffexpr*
rm -r filtered*
rm -r pathway*
#cd ../fastqc/
#rm -r *
cd ../mapped/bowtie2/featureCounts/
rm -r *
cd ../../kallisto/
rm -r *
cd ../../trimmed/
#rm -r *
cd ../
rm log.txt