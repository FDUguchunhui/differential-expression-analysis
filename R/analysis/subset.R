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


#===================================================================


#-----------------------------------------------------------------
# scenario 1
# subset a rawdata/normalized_data by a list of genes 

pos <- rownames(normalized_counts) %in% pull(gene_name)
normalized_counts_sub <- normalized_counts[pos,] 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# output is the only part that you need to give your only parameter -
#the filename
write.xlsx2(x = normalized_counts_sub, 
            file = 'output/normalized_counts_genesetsub_10hour.xlsx')

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

#--------------------------------------------------------------------------
# get full results of all comparison of interest with subset
{
  sheetname <- names(res_list_all)
  for(i in 1:length(res_list_all)){
    pos <- rownames(res_list_all[[i]]) %in% pull(gene_name)
    write.xlsx(x = res_list_all[[i]][pos,],
               file = 'output/result_10hour_genesetsub.xlsx',
               sheet = sheetname[i], append = T)
  }
}
















# get up/down regulated genes

diego <- autoHeatmap::res_subgroup(res_8hrp_2hrCF, alpha = 0.1, reg_dir = 'up')
diego2 <- autoHeatmap::res_subgroup(res_8hrp_2hrCF, alpha = 0.1, reg_dir = 'down')



#------------------------------------------------------------
# Scenario 4




upreg <- read.xlsx(file = 'data/venn_result Upregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
# fill the Names with the last Non Na value
# this step is important to get gene names of each subset
upreg <- upreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_up <- levels(factor(upreg$Names, levels = unique(upreg$Names)))

up_comb <- list()
for(i in 1:length(combo_up)){
  up_comb[[i]] <- res_all[res_all$gene_id %in% upreg[upreg$Names == combo_up[i],]$elements,]
}


#--------------------------------------------------------------------------------------------- 
for(i in 1:length(combo_up)){
  write.xlsx(up_comb[[i]],  sheetName = sub('Hours', '', combo_up)[[i]],  # the sheetname is short
             file = 'output/up_comb.xlsx', append = T)
}
#---------------------------------------------------------------------------------------------


# for down reg
downreg <- read.xlsx(file = 'data/Venn_Downregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
downreg <- downreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_down <- levels(factor(downreg$Names, levels = unique(downreg$Names)))

down_comb <- list()
for(i in 1:length(combo_down)){
  down_comb[[i]] <- res_all[res_all$gene_id %in% downreg[downreg$Names == combo_down[i],]$elements,]
}


# ------------------------------------------------------------------------------------------------
for(i in 1:length(combo_down)){
  write.xlsx(down_comb[[i]],  sheetName = gsub('hour[s]*', '', combo_down)[[i]],  # the sheetname is short
             file = 'output/down_comb.xlsx', append = T)
}
#--------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------
# scenario 5
upreg <- read.xlsx(file = 'data/venn_result Upregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
# fill the Names with the last Non Na value
# this step is important to get gene names of each subset
upreg <- upreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_up <- levels(factor(upreg$Names, levels = unique(upreg$Names)))

up_comb <- list()
for(i in 1:length(combo_up)){
  up_comb[[i]] <- normalized_counts[rownames(normalized_counts) %in% upreg[upreg$Names == combo_up[i],]$elements,]
}


#--------------------------------------------------------------------------------------------- 
for(i in 1:length(combo_up)){
  write.xlsx(up_comb[[i]],  sheetName = sub('Hours', '', combo_up)[[i]],  # the sheetname is short
             file = 'output/up_comb_normdata_kinetics.xlsx', append = T)
}
#---------------------------------------------------------------------------------------------


# for down reg
downreg <- read.xlsx(file = 'data/Venn_Downregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
# fill the Names with the last Non Na value
# this step is important to get gene names of each subset
downreg <- downreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_down <- levels(factor(downreg$Names, levels = unique(downreg$Names)))

down_comb <- list()
for(i in 1:length(combo_down)){
  down_comb[[i]] <- normalized_counts[rownames(normalized_counts) %in% downreg[downreg$Names == combo_down[i],]$elements,]
}


#--------------------------------------------------------------------------------------------- 
# sheet 7 doesn't have any observation
down_comb <- down_comb[-7]
sheetname <- gsub('hour[s]*', '', combo_down)
sheetname <- sheetname[-7]
for(i in 1:(length(combo_down) - 1)){
  write.xlsx(down_comb[[i]],  sheetName = sheetname[i],  # the sheetname is short
             file = 'output/down_comb_normdata_kinetics.xlsx', append = T)
}

