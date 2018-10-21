# SCORE

## Getting started

### Setup

Requirements:
* Anaconda
* Flexbar
* GeneSCF
* Snakemake
* Sequanix

We recommend using Anaconda since it contains many of the required libraries at installation. While it should be possible to use Miniconda, it does not contain all required packages and thus, these might need to be installed manually. Sequanix is an optional choice depending of whether the user would prefer to use a GUI or the command line.

## Usage

1.) Set the parameters in the config according to your experimental setup. If you choose to use Sequanix you can instead load the Snakefile and edit the corresponding configuration file directly from within Sequanix.

2.) Run SCORE. If you have choosen to use Sequanix you can launch the script using the GUI. If you are working on the command line execute the run-SCORE.sh script and specify the amount of cores the analysis may use. Example: ./run-SCORE.sh 8 to launch SCORE using 8 cores (default = 4).

3.) Analyze the results. All results essential for the DEG predicition are saved within the deg/ folder.

### FAQ

Q: How can I compute the runtime of my analysis?

Using the executable script mentioned above will automatically output the runtime at successful termination. Sequanix will not provide runtime information.

Q: I can't edit the config file and I can't use Sequanix, is there another way of setting parameters?

Yes, it is possible to either edit the Snakefile itself (not recommended) or using the command "snakemake -s "SCORE.snk" --config {parameter}={value}" to set {parameter} to {value}.