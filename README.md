# SCORE

## Getting started

### Setup

Requirements:
conda???
snakemake
sequanix?

## Usage

1.) Configure parameters in config directly or using sequanix
2.) Run bash script specifying core number

./SCORE.sh

to run with 4 cores (default)

./SCORE.sh <#cores>

to run with customized amount of cores

## Additional tools

### Sequanix

Using sequanix or running the snakemake file directly (without the run script) will not provide runtime information

In order to set parameters using the console, the snakefile must be run directly