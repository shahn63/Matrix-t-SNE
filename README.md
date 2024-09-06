# Matrix-t-SNE

############################################################################
##
##   t-Distributed Stochastic Neighborhood Embedding of Matrix Data
##
############################################################################

########################################
# Description of contents:
########################################

* Matrix_tSNE_v2.R

This is the main source code for Matrix t-SNE.
To run thecode, the following source codes must be loaded:

  sourceCpp("jointprob.cpp")
  sourceCpp("gradient.cpp")
  source("GaussianK.R")
  source("GaussianK_sub.R")

  
* BRCA_ex_v2.R

This script reproduces the analysis on breast cancer data.
Load the main Matrix t-SNE source code 'Matrix_tSNE.R' and the example distance matrices with:
    
  source("Matrix_tSNE_v2.R")
  load("BRCA176.RData")


* BRCA176.RData

This `RData' file contains squard distance matrices (row-wise: X1 and column-wise: X2) for 15 patients and 176 genes from breast cancer data. The gene list represents the most important genes in distinguishing BRCA1 and BRCA2, as described in Hedenfalk et al. (2001).

Reference: Hedenfalk et al. (2001). 'Gene-expression profiles in hereditary breast caner. The New Engliand Journal od Medicine.
To load the data:

  load("BRCA176.RData")


* jointprob.cpp
This C++ code in integrated with R via `Rcpp' to efficiently compute the joint probability distibution for Matrix t-SNE. By leveraging C++, the computation is significantly faster compared to pure R implementations, especially for large datasets.


* gradient.cpp
This C++ code in integrated with R via `Rcpp' to efficiently calculate the gradient of the Matrix t-SNE. Using C++ ensures that the gradient descent optimization for embedding the matrix data is performed with high computational efficiency.


* GaussianK.R
This R script provides the functionality to compute Gaussian kernels. This is a modified version of the original t-SNE code.


* GaussianJ_sub.R
This R script offers a sub-function of the Gaussian kernel computation. This is a modified versin of the original t-SNE code.


* README.txt:
  -----------

  This file provides an overview of the contents and how to use them.
  
###############################################
# R programming specifications
###############################################

It was tested with the following configuration:

    R version 4.0.5 (2021-03-31)
    Platform: x86_64-w64-mingw32/x64 (64-bit)
