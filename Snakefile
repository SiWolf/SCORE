# --------------------------------------------
# Title: SCOPE
# Author: Silver A. Wolf
# Last Modified: Thur, 14.06.2018
# Version: 0.1.2
# Usage:
#       snakemake -n
#       snakemake --dag | dot -Tsvg > dag.svg
#       snakemake -j {max_amount_of_threads}
#       snakemake --config {parameter}={value}
# --------------------------------------------

# Specifying the config file
configfile: "config.yaml"

# Global constants
PATH_BOWTIE2 = config["bowtie2_path"]
PATH_BOWTIE2_BUILD = config["bowtie2_build_path"]
PATH_FASTQC = config["fastqc_path"]
PATH_FEATURECOUNTS = config["featurecounts_path"]
PATH_FLEXBAR = config["flexbar_path"]
REF_ANNOTATION = config["ref_annotation_file"]
REF_FASTA = config["ref_fasta_file"]
REF_INDEX = config["ref_index_name"]
SAMPLES, = glob_wildcards("raw/{id}.fastq")

# Salmon
# Quasi-Mapping and transcript quanitification
# Requires a reference transcriptome as an input
# Define what is a transcriptome???
#rule salmon:
    #input:
        #"trimmed/SRR3199338_trimmed.fastq"
    #output:/home/silver/Daten/Uni/Masterarbeit 2018/Scripts/Final/Snakefile
        #"mapped/salmon/"
    #threads:
        #4
    #run:
        #shell("salmon index -t reference_transcriptome -i reference_index")
        #shell("salmon quant -i reference_index -l A -r trimmed_sample -g annotaion_file -p {threads}") #Will have to be changed for paired end sequencing    

# Master rule
# One rule to rule them all
#rule all:
    #input:
        #"deg/"

# DESeq2 Version 1.18.1
rule DEG_analysis:
    input:
        expand("mapped/bowtie2/featureCounts/{sample}/", sample=SAMPLES)
    output:
        "deg/"
    run:
        shell("Rscript libraries/SCORE.R DeSeq2 " + (" ".join(SAMPLES)))
		#shell("Rscript /libraries/SCOPE.R EdgeR SAMPLES[0] SAMPLES[1] SAMPLES[2] SAMPLES[3]")
        shell("cd ../../")
		
# featureCounts Version 1.6.2
# Counting mapped reads
# Needed for the quantification of Bowtie2 results
# Currently discards multi-mapping reads
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
# Classical genome mapping
# Always followed by quanitification/counting
rule mapping:
    input:
        "trimmed/{sample}_trimmed.fastq"
    output:
        "mapped/bowtie2/{sample}.sam"
    threads:
        4
    run:
        shell("{PATH_BOWTIE2_BUILD} {REF_FASTA} {REF_INDEX}")
        shell("{PATH_BOWTIE2} -p {threads} -x {REF_INDEX} -U {input} -S {output}")
        shell("mv {REF_INDEX}* references/")
# Verify that mapping went well?
        
# Quality control and basequality trimming
# FastQC Version 0.11.7
# Flexbar Version 3.3.0
rule quality_control_and_trimming:
    input:
        "raw/{sample}.fastq"
    output:
        "trimmed/{sample}_trimmed.fastq"
    threads:
        4
    run:
        # First FastQC run
        shell("mkdir -p fastqc/{wildcards.sample}/")
        shell("{PATH_FASTQC} {input} -o fastqc/{wildcards.sample}/")
        # Trimming
        shell("{PATH_FLEXBAR} -r {input} -t {wildcards.sample}_trimmed -n {threads} -u 20 -q TAIL -qf sanger -m 20")
        shell("mkdir -p trimmed/logs/")
        shell("mv {wildcards.sample}_trimmed.fastq trimmed/")
        shell("mv {wildcards.sample}_trimmed.log trimmed/logs/")
        shell("mkdir -p fastqc/{wildcards.sample}_trimmed/")
        # Second FastQC run
        shell("{PATH_FASTQC} trimmed/{wildcards.sample}_trimmed.fastq -o fastqc/{wildcards.sample}_trimmed/")
# Should not forget adapter trimming when using new data!