# Script for running multiple SCORE instances in different folders

cd C1/
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log.txt
cd ..
cd C2/
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log.txt
cd ..
cd C3/
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log.txt
cd ..
cd C4/
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log.txt