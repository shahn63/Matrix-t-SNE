# Matrix-t-SNE

################################################################################
##
##   t-Distributed Stochastic Neighborhood Embedding of Matrix-Variate Data
##
################################################################################

########################################
# Description of content:
########################################

* Matrix_tSNE.R

Main source code

- One needs to load the source codes 'GaussianK.R', 'GaussianK_sub.R', 'gradient.cpp', and 'jointprob.R' e.g. via

  sourceCpp("jointprob.cpp")
  sourceCpp("gradient.cpp")
  source("GaussianK.R")
  source("GaussianK_sub.R")

  
* BRCA_ex.R

Code to reproduce the analysis of the BRCA data
- One needs to load the main source code 'Matrix_tSNE.R' e.g. via
  
  source("Matrix_tSNE.R")

* README.txt:
  -----------

  This file.
  
########################################
# R programming specifications
########################################

Session info:

It was tested with the following configuration:

    R version 4.0.5 (2021-03-31)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
