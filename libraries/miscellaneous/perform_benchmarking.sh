# Script for benchmarking SCORE
# TO-DO: Split script into data generation and evaluation?
# Goal: Simulate 1x, analyse 20x
# Before: Simulate: 8x, analyse, 1x

for benchmark in {1..8}
do
	./empty_results.sh full
	cd ../../
	python libraries/fetch_transcript_lengths.py -f references/human_transduced/PROKKA_07132018.gff -i ID
	cd libraries/miscellaneous/
	
	if [ $benchmark -eq 1 ]
	then
		Rscript generate_rna_seq_data.R 50 TRUE 100 20 4
	fi
	
	if [ $benchmark -eq 2 ]
	then
		Rscript generate_rna_seq_data.R 250 TRUE 100 20 4
	fi
	
	if [ $benchmark -eq 3 ]
	then
		Rscript generate_rna_seq_data.R 50 FALSE 100 20 4
	fi
	
	if [ $benchmark -eq 4 ]
	then
		Rscript generate_rna_seq_data.R 50 TRUE 150 20 4
	fi
	
	if [ $benchmark -eq 5 ]
	then
		Rscript generate_rna_seq_data.R 50 TRUE 100 50 4
	fi
	
	if [ $benchmark -eq 6 ]
	then
		Rscript generate_rna_seq_data.R 100 FALSE 125 25 4
	fi

	if [ $benchmark -eq 7 ]
	then
		Rscript generate_rna_seq_data.R 50 TRUE 100 20 3
	fi
	
	if [ $benchmark -eq 8 ]
	then
		Rscript generate_rna_seq_data.R 50 TRUE 100 20 5
	fi

	cd simulation_data/
	source activate score_map_env
	kallisto index -i index.idx PROKKA_07132018.ffn
	for value in {1..6}
	do
		break_symbol="/"
		kallisto_name="kallisto/"
		no_1="_1.fasta.gz"
		no_2="_2.fasta.gz"
		sample="sample_0"
		sample_name_1="$sample$value$no_1"
		sample_name_2="$sample$value$no_2"
		sample_folder="$kallisto_name$sample$value$break_symbol"
		kallisto quant -i index.idx -o $sample_folder -b 100 -t 8 $sample_name_1 $sample_name_2
	done
	source deactivate
	cd ../../../
	source activate score_deg_env
	Rscript libraries/SCORE.R Metadata_C1.tsv 5000 FALSE 0.05 5 0.5 1.0 0.75 0.75 2.5 0.5 TRUE TRUE FALSE 0.5
	source deactivate
	cd libraries/miscellaneous/
	./empty_results.sh full $benchmark
done
Rscript generate_benchmarking_plots.R