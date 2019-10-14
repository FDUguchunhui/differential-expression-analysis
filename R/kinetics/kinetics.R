#this option assigns more RAM for java to enable high-RAM-comsuming processing
# especially in writing large dataset to .xlxs
options(java.parameters = "-Xmx4048m")


# # install necessary package
library(BiocManager)
# if you have install DESeq2, uncomment the following line
# BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(xlsx)
library(pheatmap)
# run next line if you need to get some help from DESeq2
# browseVignettes("DESeq2")



# read in the dataset
table <- read.csv('data/Kinetics_All_HD.csv', stringsAsFactors = F)

# 
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listFilters(mart)

# G_list_1 <- getBM(filters= "hgnc_symbol",
#                 attributes= c("ensembl_gene_id", 'hgnc_symbol')
#                 ,values=genes,mart= mart)



#--------------------------------------------------------------------------
# create the dictionary used to translate from ensemble_gene_id to 
#!!!! if it is you first use it, run it
# it will save the dataset into a R object, if you dont save R environment

#external gene_name
# G_list <- getBM(filters= "ensembl_gene_id",
#                   attributes= c("ensembl_gene_id", 'external_gene_name'),
#                   values=table$Geneid1,
#                   mart= mart)
#save(G_list, file = 'temp/G_list')

load('temp/G_list')

#-------------------------------------------------------------------
#write.xlsx2(G_list, 'gene_id_to_name.xlsx')
#------------------------------------------------------------------

#the original dataset has two gene id,
#and the second one has duplicates, in order to use the second one, need to deal with it
#!


#merge the origin with gene_name from biomart, and arrange it with original datastructure
table_merge <- sqldf::sqldf('select *
             from "table"
             left join G_list on "table".Geneid1 = G_list.ensembl_gene_id')
table_reordered <- table_merge %>%
  as_tibble() %>%
  dplyr::select(Geneid1, external_gene_name, everything())
table_reordered <- table_reordered[, c(-3, -19)]

# use the external_gene_name if exists, otherwise use ensemble_gene_id
table_reordered$external_gene_name <-
  ifelse(is.na(table_reordered$external_gene_name),table_reordered$Geneid1,
         table_reordered$external_gene_name)
# drop the Geneid1(the ensemble_id) column
table_geneid <- table_reordered[,-1]
colnames(table_geneid)[1] = 'Geneid'
# omit any row if has missing value
# Caution! this may not the best method
head(table_geneid)


#remove depublicate by geneid using SQLite
# base rule: average the read from the same gene
# caution! This may not applicable for bio-duplicate
table_nodup <-  sqldf::sqldf('
  select Geneid,
         avg(HD1_Blood) as HD1_Blood, avg(HD2_Blood) as HD2_Blood, avg(HD6_Blood) as HD6_Blood,
         avg(HD1_1H) as HD1_1H, avg(HD2_1H) as HD2_1H, avg(HD6_1H) as HD6_1H,
         avg(HD1_2H) as HD1_2H, avg(HD2_2H) as HD2_2H, avg(HD6_2H) as HD6_2H,
         avg(HD1_4H) as HD1_4H, avg(HD2_4H) as HD2_4H, avg(HD6_4H) as HD6_4H,
         avg(HD1_6H) as HD1_6H, avg(HD2_6H) as HD2_6H, avg(HD6_6H) as HD6_6H
  from table_geneid
  group by Geneid
  order by Geneid
  '
)


#check how many gene left in the processed dataset
length(table_nodup$Geneid)

#check whether there is duplicate now, if there is no duplicate, a True should be returned
all(!duplicated(table_nodup[,1]))
#check the first 5 rows
head(table_nodup)




#make the matrix needed for next step

#use the geneid as row.names 
# caution! must be unique to use as row.names(primary key)
table_rownames <- data.frame(table_nodup[,-1], row.names=table_nodup[,1])
count_maxtrix <- as.matrix(table_rownames)
# make the mode of matrix integer otherwise it will be number
mode(count_maxtrix) <- 'integer'
head(count_maxtrix)

#create the coldata corresponding to the processed dataset
coldata <- data.frame(condition = c(rep('blood', 3),
                                    rep('1h', 3),
                                    rep('2h', 3),
                                    rep('4h', 3),
                                    rep('6h', 3)),
                      row.names = colnames(count_maxtrix)[1:15])

#check whether the name of row match col, a True should be returned
all(rownames(coldata) == colnames(count_maxtrix)[1:15])

#create DESeqDataSet(dds) object, dds is a container for intermediate data
dds <- DESeqDataSetFromMatrix(countData = count_maxtrix[,1:15],
                              colData = coldata,
                              design = ~ condition)

# pre-filtering
# by removing rows in which there are very few reads, we reduce the memory size of the dds data object,
#   and we increase the speed of the transformation and testing functions within DESeq2
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set dds condition with all time points
# set 'blood' as base level, the default base level is determined by alphabet order
# other statement if also available for this purpose
# Caution!: the document only give example for factor with two levels, not sure about accuracy for
#   factor with more than 2 levels
dds$condition <- factor(dds$condition, levels = c("blood","1h",'2h', '4h', '6h'))
# drop levels that with no sample
dds$condition <- droplevels(dds$condition)



# DEseq analysis
dds <- DESeq(dds)
# for different condition you just change 6h to i.e. 1h, 2h, 4h
# may can be done by lapply() function, not sure how to set the func argument
# with different arguments
res_1h <- results(dds,contrast=c("condition","1h","blood"))
res_2h <- results(dds,contrast=c("condition","2h","blood"))
res_4h <- results(dds,contrast=c("condition","4h","blood"))
res_6h <- results(dds,contrast=c("condition","6h","blood"))

resultsNames(dds)

# log fold change shrinkage(LFC Shrinkage) may not neccessary here. for using it,
# uncomment next two lines
# resLFC <- lfcShrink(dds, coef="condition_1h_vs_blood", type="apeglm")
# resLFC


# next few lines is just a demo for using different param in res() function

# # p-values and adjusted p-values
# resOrdered <- res[order(res$pvalue),]
# summary(res)
# sum(res$padj < 0.1, na.rm=TRUE)
# res05 <- results(dds, alpha=0.05)
# summary(res05)
#
# plotMA(res, ylim=c(-2,2))
# # more useful visualize MA-plot for shrunken log2 fold changes
# plotMA(resLFC, ylim=c(-3,3))
# abline(h = c(-1,1), col='blue')
#
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]


# this function is used to subset the neeeded gene from RESULT object
#function get those LFC greate than 1 and the adjusted P-value is less than 0.1
res_subgroup <- function(res, alpha=0.1, reg_LFC=1, reg_dir='all'){

  # res is an obj from DEseq2.result function
  # alpha gives the significant level for adjusted P-value
  # reg gives the regulation level change in log2 fold change in absolute value
  # reg_dir gives which regulation direction you want to subset you gene
    # three options: all -- up and down
    #                up  -- only up regulated
    #                down -- only down regulated
  res_sig_pos <- (res$padj < alpha)
  res_sig_pos[is.na(res_sig_pos)] <- F
  if(reg_dir == 'all'){
    res_LFC_pos <- (res$log2FoldChange > reg_LFC) | (res$log2FoldChange < -reg_LFC)
  }
  else if(reg_dir == 'up'){
           res_LFC_pos <- (res$log2FoldChange > reg_LFC)
       }
       else if(reg_dir == 'down'){
                res_LFC_pos <- (res$log2FoldChange < -reg_LFC)
            }

  return(res[res_LFC_pos & res_sig_pos,])

}




