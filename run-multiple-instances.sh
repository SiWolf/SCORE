# --------------------------------
# Title: run-multiple-instances.sh
# Author: Silver A. Wolf
# Last Modified: Wed, 19.02.2020
# Version: 0.0.3
# --------------------------------

# Script for running multiple SCORE instances in a row

# Global Variables
ANALYSES=(
	I1
	I2
	I3
	I4
	M1
	M2
	M3
	M4
	W1
)
MODE="SE"
THREADS="24"

for i in "${ANALYSES[@]}"
do
	./libraries/miscellaneous/empty_results.sh full
	cp config_"$i".yaml config.yaml
	./run-SCORE.sh "$THREADS" "$MODE"
done
