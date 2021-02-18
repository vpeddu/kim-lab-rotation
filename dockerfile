# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

RUN apt-get -y update && apt-get install -y default-jdk r-cran-rjava r-cran-nloptr libssh2-1-dev

RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('viridis')"
RUN R -e "BiocManager::install('reshape2')"
RUN R -e "BiocManager::install('tidyverse')"

