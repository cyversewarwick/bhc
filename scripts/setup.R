#!/usr/bin/Rscript

update.packages()
install.packages('argparse')
install.packages('gplots')
install.packages('RColorBrewer')
source("http://bioconductor.org/biocLite.R")
biocLite("BHC")