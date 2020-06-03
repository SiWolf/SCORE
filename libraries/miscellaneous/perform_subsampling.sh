# --------------------------------
# Title: perform_benchmarking.sh
# Author: Silver A. Wolf
# Last Modified: We, 03.06.2020
# Version: 0.0.1
# --------------------------------

# Script for subsampling reads

Settings=("0.8"  "0.6"  "0.4"  "0.2")

for val in ${Settings[*]}
do
	seqtk sample -s100 raw/SRR9158217_1.fastq.gz $val > raw/$val/SRR9158217_1.fastq
	seqtk sample -s100 raw/SRR9158217_2.fastq.gz $val > raw/$val/SRR9158217_2.fastq

	seqtk sample -s100 raw/SRR9158218_1.fastq.gz $val > raw/$val/SRR9158218_1.fastq
	seqtk sample -s100 raw/SRR9158218_2.fastq.gz $val > raw/$val/SRR9158218_2.fastq

	seqtk sample -s100 raw/SRR9158219_1.fastq.gz $val > raw/$val/SRR9158219_1.fastq
	seqtk sample -s100 raw/SRR9158219_2.fastq.gz $val > raw/$val/SRR9158219_2.fastq

	seqtk sample -s100 raw/SRR9158220_1.fastq.gz $val > raw/$val/SRR9158220_1.fastq
	seqtk sample -s100 raw/SRR9158220_2.fastq.gz $val > raw/$val/SRR9158220_2.fastq

	seqtk sample -s100 raw/SRR9158222_1.fastq.gz $val > raw/$val/SRR9158222_1.fastq
	seqtk sample -s100 raw/SRR9158222_2.fastq.gz $val > raw/$val/SRR9158222_2.fastq

	seqtk sample -s100 raw/SRR9158215_1.fastq.gz $val > raw/$val/SRR9158215_1.fastq
	seqtk sample -s100 raw/SRR9158215_2.fastq.gz $val > raw/$val/SRR9158215_2.fastq
	
	gzip raw/$val/*
done