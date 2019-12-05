# SCORE: Smart Consensus Of RNA Expression pipelines

## 1.) Quick Start Guide

### 1.1) Setup

Download and extract the latest stable release of SCORE from [here](https://github.com/SiWolf/SCORE/releases).

Dependencies:
* [Anaconda](https://anaconda.org/) or [Miniconda](https://conda.io/en/latest/miniconda.html) (< 4.5.13)
* [Python 3](https://www.python.org/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (< 4.6.0)
* [tar](https://en.wikipedia.org/wiki/Tar_(computing))

Optional:
* [git](https://git-scm.com/) (Download)
* [Graphviz](https://www.graphviz.org/) (DAG visualization)
* [Sequanix](https://github.com/sequana/sequana/) (GUI)

Please ensure that the latest versions of the dependencies mentioned above are installed on your system.

Then run the installer script:

```
./libraries/miscellaneous/install.sh
```

### 1.2) Usage


1.) Feed the machine. Copy your adapter-trimmed gzipped NGS paired-end RNA-Seq reads as well as a Metadata table indicating sample groups (TSV format) into the raw/ folder. Please ensure that each read-pair shares a common prefix name (which should be consistent with the corresponding name listed in the Metadata table), followed by a suffix indicating the individual pair number (e.g. superpair_1.fastq.gz and superpair_2.fastq.gz). Please refer to [raw/Metadata.tsv](https://github.com/SiWolf/SCORE/blob/master/raw/Metadata.tsv) for an example of the required TSV format. In addition, SCORE requires a reference genome and transcriptome (both in fasta format), as well as a .gff file for annotations. These files should be placed in the references/ folder. They can be aquired for many organisms using NCBI (e.g. [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/) for E. coli) or by annotating a custom genome using tools such as [Prokka](https://github.com/tseemann/prokka). Please ensure the following steps for preprocessing your gff file: 1.) Remove the FASTA sequence located at the end of the file, 2.) Check for any unmatched " symbols and 3.) Verify that the contig names specified within your FASTA genome are matched to those within the gff file (e.g. "_size" vs. "|size").

2.) Set parameters. Edit the parameters in the configuration file ([config.yaml](https://github.com/SiWolf/SCORE/blob/master/config.yaml)) according to your experimental setup. If you choose to use a GUI, you can instead load the Snakefile using Sequanix and edit the corresponding configuration file directly from within Sequanix.

3.) Run SCORE. If you have chosen to use Sequanix you can launch SCORE using the GUI. If you are working on the command line execute the run-SCORE.sh shell script (make sure the script is executable) and specify the amount of threads the analysis may use. Example: ./run-SCORE.sh 8 to launch SCORE using 8 threads (default = 4). SCORE will automatically parallelize the internal rules accordingly.

4.) Analyze the results. All results essential for the DEG prediction are saved within the deg/<analysis_name>/ folder.

```
# Run SCORE (commandline)
# Place your raw data in raw/
# Update the config.yaml accordingly
./run-SCORE.sh <amount_of_threads>
```

For a more details please visit our [wiki](https://github.com/SiWolf/SCORE/wiki). 

## 2.) FAQ

Q: How can I compute the total runtime of my analysis?

A: Using the executable script mentioned above will automatically output the runtime at successful termination. Sequanix will not provide additional runtime information.

Q: I can't edit the config file and I can't use Sequanix, is there another way of setting parameters?

A: Yes, it is possible to either edit the Snakefile itself (not recommended) or using the command "snakemake -s "SCORE.snk" --config {parameter}={value}" to set a specified {parameter} to a {value}.

Q: Snakemake is unable to activate the included environments - what can I do?

A: This is a known issue and as of August 2019 this has not been resolved from the Snakemake developers. Following the advice mentioned [here](https://bitbucket.org/snakemake/snakemake/issues/1115/cannot-activate-conda-enironment-using), please install Conda version 4.5.13 ("conda install -n base conda=4.5.13") and modify your .bashrc (add line: export PATH="/full/path/to/miniconda/bin:$PATH") in order to ensure compability with Snakemake.

Q: I get an Snakemake error stating that not all output, log and benchmarking files contain the same wildcards? What does this mean?

A: Starting with Snakemake version 5 and above, the syntax for dynamic rules has been [updated](https://bitbucket.org/snakemake/snakemake/issues/955/problem-with-wildcard-and-dynamic-syntax). Please use Snakemake version 4.6.0 for running SCORE.

Q: How do I cite SCORE?

Please cite the following:

Wolf, S. (2018) SCORE: Smart Consensus Of RNA Expression pipelines - a consensus tool for detecting differentially expressed genes in bacteria. Free University of Berlin.