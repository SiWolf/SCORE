# Script for running multiple SCORE instances in a row

./libraries/miscellaneous/empty_results.sh full
cp config_V1.yaml config.yaml
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log_V1.txt

./libraries/miscellaneous/empty_results.sh full
cp config_V2.yaml config.yaml
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log_V2.txt

./libraries/miscellaneous/empty_results.sh full
cp config_V4.yaml config.yaml
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee log_V4.txt