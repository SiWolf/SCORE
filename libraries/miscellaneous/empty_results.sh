# Script for deleting existing SCORE results
# This is for testing various experimental combinations and parameters
# FastQC and Flexbar files should not need to be deleted

cd ../../deg/
rm -r *
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