# --------------------------------
# Title: perform_benchmarking.sh
# Author: Silver A. Wolf
# Last Modified: Thur, 12.12.2019
# Version: 0.0.7
# --------------------------------

# Script for benchmarking SCORE
# TO-DO: Split script into data generation and evaluation?
# Goal: Simulate 1x, analyze 20x
# Before: Simulate: 8x, analyze, 1x

# Global Variables
BENCHMARK_MODE="TRUE"
BOOTSTRAP="100"
INDEX="REF_INDEX"
LOW_EXPRESSION_CUTOFF="10"
MERGE_OUTPUT="TRUE"
METADATA="raw/Metadata.tsv"
NOISEQ_BIOLOGICAL_REPLICATES="FALSE"
REF_GFF="references/REF.gff"
REF_FEATURE="gene"
REF_ID="locus_tag"
REF_TRANSCRIPTOME="references/REF.ffn"
STRICT_MODE="TRUE"
SCORE_THRESHOLD="0.5"
THREADS="4"
THRESHOLD="0.05"
TOTAL_GENES="5000"
WEIGHT_BAYSEQ="1.0"
WEIGHT_DESEQ2="1.0"
WEIGHT_EDGER="1.0"
WEIGHT_LIMMA="1.0"
WEIGHT_NOISEQ="1.0"
WEIGHT_SLEUTH="1.0"

# Preprocessing
cd ../../
./libraries/miscellaneous/empty_results.sh full

# 8 Simulation Rounds
for BENCHMARK in {1..8}
do
	python3 libraries/miscellaneous/fetch_transcript_lengths.py -a $REF_GFF -f $REF_FEATURE -i $REF_ID
	python3 libraries/miscellaneous/preprocess_transcriptome.py -f $REF_TRANSCRIPTOME -i $REF_ID

	# Simulate RNA-Seq data with different parameters
	if [ $BENCHMARK -eq 1 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="50"
		RANDOMIZED="FALSE"
		READ_LENGTH="100"
	fi

	if [ $BENCHMARK -eq 2 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="100"
	fi

	if [ $BENCHMARK -eq 3 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="250"
		RANDOMIZED="TRUE"
		READ_LENGTH="100"
	fi

	if [ $BENCHMARK -eq 4 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="50"
	fi

	if [ $BENCHMARK -eq 5 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="150"
	fi

	if [ $BENCHMARK -eq 6 ]
	then
		COVERAGE="30"
		DEG_FOLD_CHANGE="4"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="100"
	fi

	if [ $BENCHMARK -eq 7 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="2"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="100"
	fi

	if [ $BENCHMARK -eq 8 ]
	then
		COVERAGE="20"
		DEG_FOLD_CHANGE="6"
		DEGS_PER_GROUP="50"
		RANDOMIZED="TRUE"
		READ_LENGTH="100"
	fi

	Rscript libraries/miscellaneous/generate_rna_seq_data.R $REF_TRANSCRIPTOME.tmp $COVERAGE $DEG_FOLD_CHANGE $DEGS_PER_GROUP $RANDOMIZED $READ_LENGTH

	# Kallisto Quantification
	source activate score_map_env
	kallisto index -i $INDEX.idx $REF_TRANSCRIPTOME.tmp

	for SAMPLE_NUMBER in {1..6}
	do
		BREAK_SYMBOL="/"
		KALLISTO_FOLDER="libraries/miscellaneous/simulation_data/kallisto/"
		PREFIX_FOLDER="libraries/miscellaneous/simulation_data/"
		RP_1="_1.fasta.gz"
		RP_2="_2.fasta.gz"
		SAMPLE="sample_0"
		SAMPLE_NAME_1="$PREFIX_FOLDER$SAMPLE$SAMPLE_NUMBER$RP_1"
		SAMPLE_NAME_2="$PREFIX_FOLDER$SAMPLE$SAMPLE_NUMBER$RP_2"
		SAMPLE_FOLDER="$KALLISTO_FOLDER$SAMPLE$SAMPLE_NUMBER$BREAK_SYMBOL"
		kallisto quant -i $INDEX.idx -o $SAMPLE_FOLDER -b $BOOTSTRAP -t $THREADS $SAMPLE_NAME_1 $SAMPLE_NAME_2
	done

	source deactivate
	rm $INDEX.idx

	# DEG Prediction
	source activate score_deg_env
	Rscript libraries/SCORE.R $METADATA $TOTAL_GENES $MERGE_OUTPUT $THRESHOLD $LOW_EXPRESSION_CUTOFF $WEIGHT_BAYSEQ $WEIGHT_DESEQ2 $WEIGHT_EDGER $WEIGHT_LIMMA $WEIGHT_NOISEQ $WEIGHT_SLEUTH $BENCHMARK_MODE $STRICT_MODE $NOISEQ_BIOLOGICAL_REPLICATES $SCORE_THRESHOLD $REF_ID
	source deactivate
	./libraries/miscellaneous/empty_results.sh full B$BENCHMARK
done

# Visualization
Rscript libraries/miscellaneous/generate_benchmarking_plots.R