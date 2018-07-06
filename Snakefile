# --------------------------------------------
# Title: SCORE
# Author: Silver A. Wolf
# Last Modified: Wed, 04.07.2018
# Version: 0.1.6
# Usage:
#		sequanix
#       snakemake -n
#       snakemake --dag | dot -Tsvg > dag.svg
#       snakemake -j {max_amount_of_threads}
#       snakemake --config {parameter}={value}
# --------------------------------------------

# Imports
import csv

# Specifying the config file
configfile: "config.yaml"

# Global parameters
METADATA = config["metadata_file"]
PATH_BOWTIE2 = config["bowtie2_path"]
PATH_BOWTIE2_BUILD = config["bowtie2_build_path"]
PATH_FASTQC = config["fastqc_path"]
PATH_FEATURECOUNTS = config["featurecounts_path"]
PATH_FLEXBAR = config["flexbar_path"]
REF_ANNOTATION = config["ref_annotation_file"]
REF_FASTA = config["ref_fasta_file"]
REF_INDEX = config["ref_index_name"]

# Functions
def read_tsv(tsv_filename):
	samples_dic = {}
	with open("raw/" + tsv_filename) as tsv:
		for line in csv.reader(tsv, delimiter = "\t"):
			if line[0][0] != "#":
				samples_dic[line[0]] = line[1]
	return(samples_dic)

# Samples
SAMPLES_AND_CONDITIONS = read_tsv(METADATA)
SAMPLES = SAMPLES_AND_CONDITIONS.keys()

#print("Welcome to SCORE: Smart Consensus Of RNA-Seq Expression pipelines")
#print("Please ensure all input files are located within the raw/ folder and any parameters have been set accordingly.")

# DESeq2 Version 1.18.1
# edgeR Version 3.20.9
rule DEG_analysis:
	input:
		expand("mapped/bowtie2/featureCounts/{sample}/", sample=SAMPLES)
	output:
		"deg/"
	run:
		# Idea: Rscript <folder>/SCORE.R <SAMPLES>
		# DESeq2 / edgeR
		shell("Rscript libraries/SCORE.R " + {METADATA})
		
# featureCounts Version 1.6.2
# Counts mapped reads to genomic features
# Needed for the quantification of Bowtie2 results
# Discards multi-mapping reads by default
rule counting:
	input:
		"mapped/bowtie2/{sample}.sam"
	output:
		"mapped/bowtie2/featureCounts/{sample}/"
	threads:
		4
	run:
		shell("{PATH_FEATURECOUNTS} -T {threads} -a {REF_ANNOTATION} -o counts {input} -g gene_name")
		shell("mv counts* {output}")

# Bowtie2 Version 2.3.4.1
# Ungapped genome mapping
# Followed by quanitification of transcripts (counting of reads)
rule mapping:
	input:
		"trimmed/{sample}_trimmed.fastq"
	output:
		"mapped/bowtie2/{sample}.sam"
	threads:
		4
	run:
		shell("{PATH_BOWTIE2_BUILD} {REF_FASTA} {REF_INDEX}")
		#shell("{PATH_BOWTIE2} -q --phred33 -p {threads} -x {REF_INDEX} -U {input} -S {output}")
		#shell("{PATH_BOWTIE2} -q --phred33 -p {threads} --no-unal -x {REF_INDEX} -1 {input_1} -2 {input_2} -S {output}")
		shell("{PATH_BOWTIE2} -q --phred33 -p {threads} --no-unal -x {REF_INDEX} -U {input} -S {output}")
		shell("mv {REF_INDEX}* references/")
# Verify that mapping went well?
        
# FastQC Version 0.11.7
# Flexbar Version 3.3.0
# Quality control and basequality trimming
rule quality_control_and_trimming:
	input:
		"raw/{sample}.fastq.gz"
	output:
		"trimmed/{sample}_trimmed.fastq"
	threads:
		4
	run:
		# Initial FastQC
		shell("mkdir -p fastqc/{wildcards.sample}/")
		shell("{PATH_FASTQC} {input} -o fastqc/{wildcards.sample}/")
		# Trimming
		shell("{PATH_FLEXBAR} -r {input} -t {wildcards.sample}_trimmed -n {threads} -u 20 -q TAIL -qf sanger -m 20")
		shell("mkdir -p trimmed/logs/")
		shell("mv {wildcards.sample}_trimmed.fastq trimmed/")
		shell("mv {wildcards.sample}_trimmed.log trimmed/logs/")
		shell("mkdir -p fastqc/{wildcards.sample}_trimmed/")
		# Second FastQC
		shell("{PATH_FASTQC} trimmed/{wildcards.sample}_trimmed.fastq -o fastqc/{wildcards.sample}_trimmed/")
# Should not forget adapter trimming when using new data!