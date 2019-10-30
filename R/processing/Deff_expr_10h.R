# don't source this script, just run-by-line first half before analysis


###! run this before you load any package. If you not, restart your session
#this option assigns more RAM for java to enable high-RAM-comsuming processing
# especially in writing large dataset to .xlxs
#memory.limit(102400)
options(java.parameters = "-Xmx10048m")

# # install necessary package
library(BiocManager)
# if you have install DESeq2, uncomment the followiing line
# BiocManager::install("DESeq2")
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(xlsx)
library(tidyverse)
library(magrittr)
library(dplyr)
# run next line if you need to get some help from DESeq2
# browseVignettes("DESeq2")




gene_expr_10h.table <- read.xlsx2('data/10hr_normed_counts.xlsx', sheetName = 'counts', stringsAsFactors = F,
                                  colClasses = c('character', rep('numeric', 16)))
gene_expr_10h.tibble <- as_tibble(gene_expr_10h.table)
#Check all unique
all(!duplicated(gene_expr_10h.tibble[,1]))


#use the geneid as row.names
table_rownames <- data.frame(gene_expr_10h.tibble[,-1], row.names=gene_expr_10h.table[,1])
count_maxtrix <- as.matrix(table_rownames)
count_maxtrix <- round(count_maxtrix)
# make the mode of matrix integer otherwise it will be number
head(count_maxtrix)

coldata <- data.frame(condition = c(rep('blood', 4),
                                    rep('LTB4_TM', 5),
                                    rep('CFASN_TM', 5),
                                    rep('CFASN_inc', 2)),
                      row.names = colnames(count_maxtrix[,1:16]))

all(rownames(coldata) == colnames(count_maxtrix[1:16]))

#create DESeqDataSet(dds) object, dds is a container for intermediate data
dds <- DESeqDataSetFromMatrix(countData = count_maxtrix[,1:16],
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
dds$condition <- factor(dds$condition, levels = c("blood","LTB4_TM",'CFASN_TM', 'CFASN_inc'))
# drop levels that with no sample
dds$condition <- droplevels(dds$condition)


# DEseq analysis
dds <- DESeq(dds)


# for different condition you just change 6h to i.e. 1h, 2h, 4h
# may can be done by lapply() function, not sure how to set the func argument
# with different arguments
res_LTB4 <- results(dds, contrast=c("condition","LTB4_TM", "blood"))
res_CFASN <- results(dds, contrast=c("condition","CFASN_TM", "blood"))

# comparison within inc and blood
res_inc <- results(dds, contrast=c("condition","CFASN_inc", "blood"))

res_list_all <- list(res_LTB4, res_CFASN, res_inc)









# do not run the following code
#--------------------------------------------------------------------------------
#upregulated
res_LTB4_up <- as.data.frame(res_subgroup(res_LTB4, reg_dir = 'up'))
res_LTB4_up$gene_name <- row.names(res_LTB4_up)
res_CFASN_up <- as.data.frame(res_subgroup(res_CFASN, reg_dir = 'up'))
res_CFASN_up$gene_name <- row.names(res_CFASN_up)
res_inc_up <- as.data.frame(res_subgroup(res_inc, reg_dir = 'up'))
res_inc_up$gene_name <- row.names(res_inc_up)

#down regulated
res_LTB4_down <- as.data.frame(res_subgroup(res_LTB4, reg_dir = 'down'))
res_LTB4_down$gene_name <- row.names(res_LTB4_down)
res_CFASN_down <- as.data.frame(res_subgroup(res_CFASN, reg_dir = 'down'))
res_CFASN_down$gene_name <- row.names(res_CFASN_down)
res_inc_down <- as.data.frame(res_subgroup(res_inc, reg_dir = 'down'))
res_inc_down$gene_name <- row.names(res_inc_down)

#write significant up & down regulated genes
write.xlsx2(res_LTB4_up,file = 'res_LTB4_up.xlsx')
write.xlsx2(res_CFASN_up, file = 'res_CFASN_up.xlsx')
write.xlsx2(res_LTB4_down,file = 'res_LTB4_down.xlsx')
write.xlsx2(res_CFASN_down, file = 'res_CFASN_down.xlsx')
write.xlsx2(res_inc_up, file = 'output/res_inc_up.xlsx')
write.xlsx2(res_inc_down, file = 'output/res_inc_down.xlsx')

# intersection for up and down
res_up_intersection <- sqldf::sqldf('
              select *
              from res_LTB4_up
              inner join res_CFASN_up on
                  res_LTB4_up.gene_name = res_CFASN_up.gene_name
             ') %>%   dplyr::select(gene_name, baseMean:padj)

res_down_intersection <- sqldf::sqldf('
              select *
              from res_LTB4_down
              inner join res_CFASN_down on
                  res_LTB4_down.gene_name = res_CFASN_down.gene_name
             ') %>%   dplyr::select(gene_name, baseMean:padj)

write.xlsx2(res_up_intersection,file = 'res_up_intersection.xlsx')
write.xlsx2(res_down_intersection,file = 'res_down_intersection.xlsx')





# take unique
res_LTB4_up.unique <- sqldf::sqldf('
              select *
              from res_LTB4_up
              where res_LTB4_up.gene_name not in (
                select res_LTB4_up.gene_name
                from res_LTB4_up
                inner join res_CFASN_up on
                  res_LTB4_up.gene_name = res_CFASN_up.gene_name
              )
             ')

res_CFASN_up.unique <- sqldf::sqldf('
              select *
              from res_CFASN_up
              where res_CFASN_up.gene_name not in (
                select res_CFASN_up.gene_name
                from res_CFASN_up
                inner join res_LTB4_up on
                  res_LTB4_up.gene_name = res_CFASN_up.gene_name
              )
             ')

#dwon regulated
res_LTB4_down.unique <- sqldf::sqldf('
              select *
              from res_LTB4_down
              where res_LTB4_down.gene_name not in (
                select res_LTB4_down.gene_name
                from res_LTB4_down
                inner join res_CFASN_down on
                  res_LTB4_down.gene_name = res_CFASN_down.gene_name
              )
             ')

res_CFASN_down.unique <- sqldf::sqldf('
              select *
              from res_CFASN_down
              where res_CFASN_down.gene_name not in (
                select res_CFASN_down.gene_name
                from res_CFASN_down
                inner join res_LTB4_down on
                  res_LTB4_down.gene_name = res_CFASN_down.gene_name
              )
             ')

res_LTB4_up.unique %<>% as_tibble()
res_LTB4_up.unique <- res_LTB4_up.unique[,1:7] %>% dplyr::select(gene_name, everything())
res_CFASN_up.unique <- res_CFASN_up.unique[,1:7] %>% dplyr::select(gene_name, everything())

#down
res_LTB4_down.unique <- res_LTB4_down.unique[,1:7] %>% select(gene_name, everything())
res_CFASN_down.unique <- res_CFASN_down.unique[,1:7] %>% select(gene_name, everything())

#up
write.xlsx2(res_LTB4_up.unique,file = 'res_LTB4_up_unique.xlsx')
write.xlsx2(res_CFASN_up.unique,file = 'res_CFASN_up_unique.xlsx')

#down regulated
write.xlsx2(res_LTB4_down.unique,file = 'res_LTB4_down_unique.xlsx')
write.xlsx2(res_CFASN_down.unique,file = 'res_CFASN_down_unique.xlsx')



# check unique
all(!(res_LTB4_up.unique$gene_name %in% res_CFASN_up.unique$gene_name))
all(!(res_LTB4_down.unique$gene_name %in% res_CFASN_down.unique$gene_name))

# unique 10h up-regulated with other 4 time points CFASN
length(res_CFASN_up.unique$gene_name)
intersc_unique_10_CFASN_up_four_timepoints <- res_CFASN_up.unique[res_CFASN_up.unique$gene_name %in% genes_up_intersection,]
length(intersc_unique_10_CFASN_up_four_timepoints$gene_name)
write.xlsx2(intersc_unique_10_CFASN_up_four_timepoints, 'intersc_unique_10_CFASN_up_four_timepoints.xlsx')

#make a table with the lfc of 5 time points of up-regulated CFASN
# intersc_unique_10_CFASN_up_four_timepoints.lfc_ts_expr.tibble <- lfc_ts_expr.tibble %>%
#   filter(gene_name %in% intersc_unique_10_CFASN_up_four_timepoints$gene_name)
#
# intersc_unique_10_CFASN_up_four_timepoints.lfc_ts_expr.tibble <- sqldf::sqldf('
#   select
#     "lfc.1h", "lfc.2h", "lfc.4h", "lfc.6h", log2FoldChange as "lfc.10h"
#     from "intersc_unique_10_CFASN_up_four_timepoints.lfc_ts_expr.tibble", "res_CFASN_up.unique"
#     where "intersc_unique_10_CFASN_up_four_timepoints.lfc_ts_expr.tibble".gene_name = "res_CFASN_up.unique".gene_name
# ')
#
#
# write.table(intersc_unique_10_CFASN_up_four_timepoints.lfc_ts_expr.tibble, file = 'intersc_unique_10_CFASN_up_four_timepoints.tsv',
#             quote = F, sep = '\t', row.names = F)

intersc_unique_10_CFASN_up_four_timepoints


