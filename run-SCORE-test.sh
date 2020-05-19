# ------------------------------
# Title: run-SCORE-test.sh
# Author: Silver A. Wolf
# Last Modified: Wed, 19.05.2020
# Version: 0.0.1
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
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
	gunzip references/*
	snakemake -s "SCORE_PE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
	{ time snakemake -s "SCORE_PE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
fi

if [ $STUDY -eq 2 ]
then
	mv config_Peyrusson_et_al.yaml config.yaml
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz
	snakemake -s "SCORE_PE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
	{ time snakemake -s "SCORE_PE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
fi

if [ $STUDY -eq 3 ]
then
	mv config_Rodman_et_al.yaml config.yaml
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.fna.gz
	wget -P references/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/746/645/GCF_000746645.1_ASM74664v1/GCF_000746645.1_ASM74664v1_genomic.gff.gz
	gunzip references/*
	snakemake -s "SCORE_PE.snk" --forceall --dag | dot -Tsvg > deg/dag.svg
	{ time snakemake -s "SCORE_PE.snk" -j $NUM_CORES --use-conda ; } 2>&1 | tee deg/log.txt
fi

mv config.yaml.tmp config.yaml