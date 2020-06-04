
# This container will install SCORE from master
#
FROM continuumio/miniconda3
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda install -y python=3.6 snakemake=4.6.0 gffutils



#Get the latest code from github and install
RUN git clone https://github.com/SiWolf/SCORE
RUN ["chmod", "+x", "/SCORE/libraries/miscellaneous/install.sh"]
RUN ["chmod", "+x", "/SCORE/run-SCORE.sh"]

RUN cd /SCORE && ./libraries/miscellaneous/install.sh
