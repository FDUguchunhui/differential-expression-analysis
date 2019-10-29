#this option assigns more RAM for java to enable high-RAM-comsuming processing
# especially in writing large dataset to .xlxs
options(java.parameters = "-Xmx8048m")


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
library(xlsx)

# get the target list first
#==================================================================
## subset gene
gene_name <- read_table(file = 'data/geneset.txt', skip = 2, 
                        col_names = 'gene name',
                        col_types = cols(
                          `gene name` = col_character()
                        ))

pos <- rownames(normalized_counts) %in% pull(gene_name)
#===================================================================


#-----------------------------------------------------------------
# scenario 1
# subset a rawdata/normalized_data by a list of genes 


normalized_counts_sub <- normalized_counts[pos,] 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# output is the only part that you need to give your only parameter -
#the filename
write.xlsx2(x = normalized_counts_sub, 
            file = 'output/normalized_counts_sub.xlsx')

#------------------------------------------------------------------------


#------------------------------------------------------------------------
# scenario 2
# result_from_DESeq2




# kinetics_all
{
  res_list <-
    list(
      res_1h = res_1h,
      res_2h = res_2h,
      res_4h = res_4h,
      res_6h = res_6h
    )
  DEresult(res_list)

}

#--------------------------------------------------------------------
# get up/down regulated genes

diego <- autoHeatmap::res_subgroup(res_8hrp_2hrCF, alpha = 0.1, reg_dir = 'up')
diego2 <- autoHeatmap::res_subgroup(res_8hrp_2hrCF, alpha = 0.1, reg_dir = 'down')



#------------------------------------------------------------
