# SCORE: Smart Consensus Of RNA Expression pipelines

## 1.) Getting started

### 1.1) Setup

Download and extract the latest stable release of SCORE from [here](https://github.com/SiWolf/SCORE/releases).

Dependencies:
* [Anaconda](https://anaconda.org/)
* [Flexbar](https://github.com/seqan/flexbar)
* [GeneSCF](http://genescf.kandurilab.org/)
* [Python3](https://www.python.org/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)

Optional:
* [Sequanix](https://github.com/sequana/sequana/) (Part of Sequana)

Please ensure that the latest versions of the dependencies mentioned above are installed on your system. For Conda environments, it is recommend to use Anaconda since it contains many of the required libraries at installation. While it should be possible to use Miniconda, it does not contain all required packages and thus, these might need to be installed manually. Sequanix is an optional choice depending on whether or not the user would prefer to use a GUI or prefers to work on the command line. Sequanix should not be installed for the latter case. SCORE primarily supports Ubuntu Linux distributions but might work on other systems as long as the dependencies mentioned above are installed accordingly.

### 1.2) Usage

1.) Feed the machine. Copy your adapter-trimmed NGS paired-end RNA-Seq reads as well as a Metadata table indicating sample groups (TSV format) into the raw/ folder. Please ensure that each read-pair shares a common prefix name (which should be consistent with the corresponding name listed in the Metadata table), followed by a suffix indicating the individual pair number (e.g. superpair_1.fastq.gz and superpair_2.fastq.gz). Please refer to raw/Metadata_old.tsv for an example of the required TSV format.

2.) Set parameters. Edit the parameters in the configuration file (config.yaml) according to your experimental setup. If you choose to use a GUI, you can instead load the Snakefile using Sequanix and edit the corresponding configuration file directly from within Sequanix.

3.) Run SCORE. If you have chosen to use Sequanix you can launch SCORE using the GUI. If you are working on the command line execute the run-SCORE.sh shell script (make sure the script is executable) and specify the amount of threads the analysis may use. Example: ./run-SCORE.sh 8 to launch SCORE using 8 threads (default = 4). SCORE will automatically parallelize the internal rules accordingly.

4.) Analyze the results. All results essential for the DEG prediction are saved within the deg/<analysis_name>/ folder.

## 2.) FAQ

Q: How can I compute the total runtime of my analysis?

Using the executable script mentioned above will automatically output the runtime at successful termination. Sequanix will not provide additional runtime information.

Q: I can't edit the config file and I can't use Sequanix, is there another way of setting parameters?

Yes, it is possible to either edit the Snakefile itself (not recommended) or using the command "snakemake -s "SCORE.snk" --config {parameter}={value}" to set a specified {parameter} to a {value}.

Q: How do I cite SCORE?

Please cite the following:

Wolf, S. (2018) SCORE: Smart Consensus Of RNA Expression pipelines - a consensus tool for detecting differentially expressed genes in bacteria. Free University of Berlin.