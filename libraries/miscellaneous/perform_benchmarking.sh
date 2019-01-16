# Script for benchmarking SCORE

for benchmark in {1..5}
do
	./empty_results.sh full
	cd ../../
	python libraries/fetch_transcript_lengths.py -f references/human_transduced/PROKKA_07132018.gff -i ID
	cd libraries/miscellaneous/
	Rscript generate_rna_seq_data.R
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
	cd ../../
	source activate score_deg_env
	Rscript SCORE.R
	source deactivate
	cd miscellaneous/
	./empty_results.sh full $benchmark
done