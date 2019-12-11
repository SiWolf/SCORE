# SCORE
>Smart Consensus Of RNA Expression pipelines

### Table of Contents

1. [Setup](#1-setup)
2. [Usage](#2-usage)
3. [FAQ](#3-faq)
4. [License](#4-license)

## 1. Setup

Download and extract the latest stable release of SCORE from [here](https://github.com/SiWolf/SCORE/releases).

**Dependencies:**
* [Anaconda](https://anaconda.org/) or [Miniconda](https://conda.io/en/latest/miniconda.html) (< 4.5.13)
* [Python 3](https://www.python.org/)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/) (< 4.6.0)
* [tar](https://en.wikipedia.org/wiki/Tar_(computing))

**Optional:**
* [git](https://git-scm.com/) (Download)
* [Graphviz](https://www.graphviz.org/) (DAG visualization)
* [Sequanix](https://github.com/sequana/sequana/) (GUI)

Please ensure that the latest versions of the dependencies mentioned above are installed on your system.

Then run the installation script:

```
./libraries/miscellaneous/install.sh
```

## 2. Usage

**Requirements:**
* Adapter-Trimmed NGS Reads (.fastq.gz)
* Reference Annotation File (.gff)
* Reference Genome (.fasta)
* Reference Transcriptome (.fasta)

Place input reads into the [raw/](https://github.com/SiWolf/SCORE/tree/master/raw) directory and the reference files into [references/](https://github.com/SiWolf/SCORE/tree/master/references).

Set parameters in the following files:
* [config.yaml](https://github.com/SiWolf/SCORE/blob/master/config.yaml)
* [raw/Metadata.tsv](https://github.com/SiWolf/SCORE/blob/master/raw/Metadata.tsv)

And run SCORE:

```
./run-SCORE.sh <amount_of_threads>
```

For more details please refer to our [wiki](https://github.com/SiWolf/SCORE/wiki).

## 3. FAQ

**How can I compute the total runtime of my analysis?**

Using the executable script mentioned above will automatically output the runtime at successful termination. Sequanix will not provide additional runtime information.

**I can't edit the config file and I can't use Sequanix, is there another way of setting parameters?**

Yes, it is possible to either edit the Snakefile itself (not recommended) or using the command "snakemake -s "SCORE.snk" --config {parameter}={value}" to set a specified {parameter} to a {value}.

**Snakemake is unable to activate the included environments - what can I do?**

This is a known issue and as of August 2019 this has not been resolved from the Snakemake developers. Following the advice mentioned [here](https://bitbucket.org/snakemake/snakemake/issues/1115/cannot-activate-conda-enironment-using), please install Conda version 4.5.13 ("conda install -n base conda=4.5.13") and modify your .bashrc (add line: export PATH="/full/path/to/miniconda/bin:$PATH") in order to ensure compability with Snakemake.

**Snakemake outputs an error stating that not all output, log and benchmarking files contain the same wildcards? What does this mean?**

Starting with Snakemake version 5 and above, the syntax for dynamic rules has been [updated](https://bitbucket.org/snakemake/snakemake/issues/955/problem-with-wildcard-and-dynamic-syntax). Please use Snakemake version 4.6.0 for running SCORE.

**How do I cite SCORE?**

Please cite the following:

Wolf, S. (2018) SCORE: Smart Consensus Of RNA Expression pipelines - a consensus tool for detecting differentially expressed genes in bacteria. Free University of Berlin.

## 4. License

This project is licensed under the GPLv3 License. See the [LICENSE.md](LICENSE) file for details.
