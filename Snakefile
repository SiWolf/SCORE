# --------------------------------------------
# Title: SCORE
# Author: Silver A. Wolf
# Last Modified: Thur, 09.08.2018
# Version: 0.2.2
# Usage:
#		snakemake -j {max_amount_of_threads} snakemake --use-conda
# Additional options:
#		sequanix
#       snakemake -n
#       snakemake --dag | dot -Tsvg > dag.svg
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
PATH_FLEXBAR = config["flexbar_path"]
REF_ANNOTATION = config["ref_annotation_file"]
REF_ANNOTATION_FEATURE_ID = config["ref_annotation_feature_type"]
REF_ANNOTATION_GENE_ID = config["ref_annotation_gene_type"]
REF_FASTA = config["ref_fasta_file"]
REF_INDEX = config["ref_index_name"]

# Functions
def read_tsv(tsv_filename):
	samples_dic = {}
	with open("raw/" + tsv_filename) as tsv:
		for line in csv.reader(tsv, delimiter = "\t"):
			if line[0][0] != "@":
				samples_dic[line[0]] = line[1]
	return(samples_dic)

# Samples
SAMPLES_AND_CONDITIONS = read_tsv(METADATA)
SAMPLES = SAMPLES_AND_CONDITIONS.keys()

rule postprocessing:
	input:
		"deg_analysis_graphs.pdf"
	output:
		"deg/deg_analysis_graphs.pdf"
	run:
		# Moving the alignment reference index file now, since it's not used anymore
		shell("mv {REF_INDEX}* references/")
		shell("mv {input} deg/")

# baySeq Version 2.12.0
# DESeq2 Version 1.18.1
# edgeR Version 3.20.7
rule DEG_analysis:
	input:
		expand("mapped/bowtie2/featureCounts/{sample}/", sample = SAMPLES)
	output:
		"deg_analysis_graphs.pdf"
	conda:
		"libraries/score_deg_environment.yml"
	shell:
		# Idea: Rscript <folder>/SCORE.R <SAMPLES>
		"Rscript libraries/SCORE.R {METADATA}"
		
# featureCounts Version 1.6.2
# Counts mapped reads to genomic features
# Needed for the quantification of Bowtie2 results
# Discards multi-mapping reads by default
rule counting:
	input:
		"mapped/bowtie2/{sample}.sam"
	output:
		"mapped/bowtie2/featureCounts/{sample}/"
	conda:
		"libraries/score_count_environment.yml"
	threads:
		4
	shell:
		"featureCounts -T {threads} -a {REF_ANNOTATION} -o counts {input} -t {REF_ANNOTATION_FEATURE_ID} -g {REF_ANNOTATION_GENE_ID} && mv counts* {output}"

# Bowtie2 Version 2.3.4.1
# Ungapped genome mapping
# Followed by quanitification of transcripts (counting of reads)
# TO-DO: Verify that mapping went well?
# TO-DO: Create Conda environment for Bowtie2 (does not work as of 9-8-18)
rule mapping:
	input:
		"trimmed/{sample}_trimmed_1.fastq.gz",
		"trimmed/{sample}_trimmed_2.fastq.gz"
	output:
		"mapped/bowtie2/{sample}.sam"
	threads:
		4
	run:
		shell("{PATH_BOWTIE2_BUILD} {REF_FASTA} {REF_INDEX}")
		p1 = input[0]
		p2 = input[1]
		shell("{PATH_BOWTIE2} -q --phred33 -p {threads} --no-unal -x {REF_INDEX} -1 {p1} -2 {p2} -S {output}")
        
# FastQC Version 0.11.7
# Flexbar Version 3.3.0
# Quality control and basequality trimming
# TO-DO: Should not forget adapter trimming when using new data!
# TO-DO: Create Conda environment for FastQC and Flexbar (does not work as of 9-8-18)
rule quality_control_and_trimming:
	input:
		"raw/{sample}_1.fastq.gz",
		"raw/{sample}_2.fastq.gz"
	output:
		"trimmed/{sample}_trimmed_1.fastq.gz",
		"trimmed/{sample}_trimmed_2.fastq.gz"
	log:
		"trimmed/logs/"
	threads:
		4
	run:
		# Initial FastQC
		shell("mkdir -p fastqc/{wildcards.sample}/")
		shell("{PATH_FASTQC} {input} -o fastqc/{wildcards.sample}/")
		# Trimming
		r1 = input[0]
		r2 = input[1]
		shell("{PATH_FLEXBAR} -r {r1} -p {r2} -t {wildcards.sample}_trimmed -n {threads} -u 20 -q TAIL -qf sanger -m 20 -z GZ")
		shell("mv {wildcards.sample}_trimmed_1.fastq.gz trimmed/")
		shell("mv {wildcards.sample}_trimmed_2.fastq.gz trimmed/")
		shell("mv {wildcards.sample}_trimmed.log {log}")
		shell("mkdir -p fastqc/{wildcards.sample}_trimmed/")
		# Second FastQC
		shell("{PATH_FASTQC} trimmed/{wildcards.sample}_trimmed_1.fastq.gz -o fastqc/{wildcards.sample}_trimmed/")
		shell("{PATH_FASTQC} trimmed/{wildcards.sample}_trimmed_2.fastq.gz -o fastqc/{wildcards.sample}_trimmed/")

onstart:
	print("\n Welcome to SCORE: Smart Consensus Of RNA-Seq Expression pipelines")
	print(" Please ensure all input files are located within the raw/ folder and your parameters have been set accordingly.")		
		
onerror:
	print("\n Something went wrong. Please refer to the error messages above. \n")

onsuccess:
	print("\n SCORE finished successfully. \n")