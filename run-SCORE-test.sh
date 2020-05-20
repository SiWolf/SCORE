# ------------------------------
# Title: run-SCORE-test.sh
# Author: Silver A. Wolf
# Last Modified: Wed, 20.05.2020
# Version: 0.0.3
# ------------------------------

# Script for performing SCORE for 3 test datasets
# Selection of different studies:
# 	(1) Houser et. al. 
# 	(2) Peyrusson et. al.
# 	(3) Rodman et. al.
# How to run: ./run-SCORE.sh <#cores> <study>

NUM_CORES=$1
STUDY=$2

./libraries/miscellaneous/empty_results.sh full

mv config.yaml config.yaml.tmp

if [ $STUDY -eq 1 ]
then
	mv config_Houser_et_al.yaml config.yaml
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/001/SRR1945181/SRR1945181_1.fastq.gz -o raw/SRR1945181_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/001/SRR1945181/SRR1945181_2.fastq.gz -o raw/SRR1945181_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/006/SRR1945196/SRR1945196_1.fastq.gz -o raw/SRR1945196_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/006/SRR1945196/SRR1945196_2.fastq.gz -o raw/SRR1945196_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/007/SRR1945197/SRR1945197_1.fastq.gz -o raw/SRR1945197_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/007/SRR1945197/SRR1945197_2.fastq.gz -o raw/SRR1945197_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/006/SRR1945206/SRR1945206_1.fastq.gz -o raw/SRR1945206_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/006/SRR1945206/SRR1945206_2.fastq.gz -o raw/SRR1945206_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/009/SRR1945209/SRR1945209_1.fastq.gz -o raw/SRR1945209_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/009/SRR1945209/SRR1945209_2.fastq.gz -o raw/SRR1945209_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/004/SRR1945234/SRR1945234_1.fastq.gz -o raw/SRR1945234_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/004/SRR1945234/SRR1945234_2.fastq.gz -o raw/SRR1945234_2.fastq.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
	gunzip references/*
	{ time snakemake -s "SCORE_PE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
	mv config.yaml config_Houser_et_al.yaml
fi

if [ $STUDY -eq 2 ]
then
	mv config_Peyrusson_et_al.yaml config.yaml
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/021/SRR10379721/SRR10379721.fastq.gz -o raw/SRR10379721.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/022/SRR10379722/SRR10379722.fastq.gz -o raw/SRR10379722.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/023/SRR10379723/SRR10379723.fastq.gz -o raw/SRR10379723.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/024/SRR10379724/SRR10379724.fastq.gz -o raw/SRR10379724.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/025/SRR10379725/SRR10379725.fastq.gz -o raw/SRR10379725.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/026/SRR10379726/SRR10379726.fastq.gz -o raw/SRR10379726.fastq.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz
	gunzip references/*
	{ time snakemake -s "SCORE_SE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
	mv config.yaml config_Peyrusson_et_al.yaml
fi

if [ $STUDY -eq 3 ]
then
	mv config_Rodman_et_al.yaml config.yaml
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/007/SRR9158217/SRR9158217_1.fastq.gz -o raw/SRR9158217_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/007/SRR9158217/SRR9158217_2.fastq.gz -o raw/SRR9158217_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/008/SRR9158218/SRR9158218_1.fastq.gz -o raw/SRR9158218_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/008/SRR9158218/SRR9158218_2.fastq.gz -o raw/SRR9158218_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/009/SRR9158219/SRR9158219_1.fastq.gz -o raw/SRR9158219_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/009/SRR9158219/SRR9158219_2.fastq.gz -o raw/SRR9158219_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/000/SRR9158220/SRR9158220_1.fastq.gz -o raw/SRR9158220_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/000/SRR9158220/SRR9158220_2.fastq.gz -o raw/SRR9158220_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/001/SRR9158221/SRR9158221_1.fastq.gz -o raw/SRR9158221_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/001/SRR9158221/SRR9158221_2.fastq.gz -o raw/SRR9158221_2.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/005/SRR9158215/SRR9158215_1.fastq.gz -o raw/SRR9158215_1.fastq.gz
	curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR915/005/SRR9158215/SRR9158215_2.fastq.gz -o raw/SRR9158215_2.fastq.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.gff.gz
	gunzip references/*
	{ time snakemake -s "SCORE_PE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
	mv config.yaml config_Rodman_et_al.yaml
fi

mv config.yaml.tmp config.yaml