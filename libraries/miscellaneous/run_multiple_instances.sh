# --------------------------------
# Title: run_multiple_instances.sh
# Author: Silver A. Wolf
# Last Modified: Thue, 18.02.2020
# Version: 0.0.2
# --------------------------------

# Script for running multiple SCORE instances in a row

# Global Variables
MODE="SE"
THREADS="12"

./libraries/miscellaneous/empty_results.sh full
cp config_R1.yaml config.yaml
snakemake -s "SCORE_"$MODE".snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE_"$MODE".snk" -j $THREADS --use-conda ; } 2>&1 | tee deg/log.txt

./libraries/miscellaneous/empty_results.sh full
cp config_R2.yaml config.yaml
snakemake -s "SCORE_"$MODE".snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE_"$MODE".snk" -j $THREADS --use-conda ; } 2>&1 | tee deg/log.txt

./libraries/miscellaneous/empty_results.sh full
cp config_R3.yaml config.yaml
snakemake -s "SCORE_"$MODE".snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE_"$MODE".snk" -j $THREADS --use-conda ; } 2>&1 | tee deg/log.txt

./libraries/miscellaneous/empty_results.sh full
cp config_R4.yaml config.yaml
snakemake -s "SCORE_"$MODE".snk" --forceall --dag | dot -Tsvg > deg/dag.svg
{ time snakemake -s "SCORE_"$MODE".snk" -j $THREADS --use-conda ; } 2>&1 | tee deg/log.txt
