# Matrix-t-SNE

################################################################################
##
##   t-Distributed Stochastic Neighborhood Embedding of Matrix-Variate Data
##
################################################################################

########################################
# Description of content:
########################################

* Matrix_tSNE_v2.R

Main source code

- One needs to load the source codes 'GaussianK.R', 'GaussianK_sub.R', 'gradient.cpp', and 'jointprob.R' e.g. via

  sourceCpp("jointprob.cpp")
  sourceCpp("gradient.cpp")
  source("GaussianK.R")
  source("GaussianK_sub.R")

  
* BRCA_ex_v2.R

Code to reproduce the analysis of the BRCA data
- One needs to load the main source code 'Matrix_tSNE.R' e.g. via
  
  source("Matrix_tSNE_v2.R")

* BRCA176.RData

RData file containing row-wise (X1) and column-wise (X2) squared distance matrices of 15 patients and 176 genes from breast cancer data.
The most importance 176 genes; list in distinguishing BRCA1 and BRCA2 are available in Hedenfalk et al. (2001)
Reference: Hedenfalk et al. (2001). 'Gene-expression profiles in hereditary breast caner. The New Engliand Journal od Medicine.
- One needs to load the RData file 'BRCA176.RData' e.g. via

  load("BRCA176.RData")


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
