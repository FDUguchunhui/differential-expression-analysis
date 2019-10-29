# this script is design to used with codes in processing file
# please run kinetics.R or equivalence before running
# this script is not designed to run automatically, just for personal, quick and dirty analysis
# My suggestion is that you should check the code, then modify it or write your own

# this analysis script can be reused by many DE analysis result
#!!! my variable names are reusable, my suggestion is to clean your environment 
# and rebuild you environment before you run analysis
library(tidyverse)
# autoHeatmap is a in-house package, a sourced package is available 
#in the root directory
library(autoHeatmap)

# there are two ways to get normalized data from DESeq2

#--------------------------------------------------------------------
# method 1, if you just want the normalized data,
#you don't have to run DESeq() function
# just create normalized count
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#--------------------------------------------------------------------



#--------------------------------------------------------------------
# method 2, if you already have result from DESeq(), you can use
#following code
normalized_counts <- counts(dds, normalized=TRUE)
mean(normalized_counts[,"HD1_8hr_CFplus"] >  normalized_counts[,"HD1_8hr_CFminus"])
#---------------------------------------------------------------------



#
