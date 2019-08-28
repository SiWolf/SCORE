# Script for running multiple SCORE instances in a row

./libraries/miscellaneous/empty_results.sh full
cp config_V1.yaml config.yaml
snakemake -s "SCORE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee deg/log.txt

./libraries/miscellaneous/empty_results.sh full
cp config_V2.yaml config.yaml
snakemake -s "SCORE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee deg/log.txt

./libraries/miscellaneous/empty_results.sh full
cp config_V4.yaml config.yaml
snakemake -s "SCORE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE.snk" -j 12 --use-conda ; } 2>&1 | tee deg/log.txt