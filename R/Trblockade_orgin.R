options(java.parameters = "-Xmx10048m")

library(DESeq2)
library(tidyverse)
library(xlsx)
library(magrittr)

#---------------------------------------------------------------------
dat <- read_csv(file = 'data/Trblockade_orgin.csv',
                col_types = cols(
                  `GENE ID` = col_character(),
                  HD1_2hr_CF = col_double(),
                  HD1_8hr_CFminus = col_double(),
                  HD1_8hr_CFplus = col_double(),
                  HD2_2hr_CF = col_double(),
                  HD2_2hr_LTB4 = col_double(),
                  HD2_8hr_CFminus = col_double(),
                  HD3_2hr_CF = col_double(),
                  HD3_2hr_LTB4 = col_double(),
                  HD3_8hr_CFminus = col_double(),
                  HD3_8hr_CFplus = col_double()
                ))

# read the blank rows
dat %<>% na.omit()

dat <- dat %>% group_by(`GENE ID`) %>% 
  summarise_all(.funs = mean)

# check unique
all(!duplicated(dat$`GENE ID`))

#------------------------------------------------------------------------------
# step: create coldata
colnames(dat)
coldata <- data.frame(condition = sub('HD\\d{1}_', '', colnames(dat[-1])),
                  row.names = colnames(dat[-1]))

# create a integer matrix of count
count_matrix <- as.matrix(dat[,-1])
rownames(count_matrix) <- pull(dat[,1])
mode(count_matrix) <- 'integer'

# check whether condition in coldata is the same as count_matrix 
all(rownames(coldata) == colnames(count_matrix))



# -----------------------------------------------------------------------------
#create DESeqDataSet(dds) object, dds is a container for intermediate data
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
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
#run the one of following line
dds$condition <- factor(dds$condition, levels = c("2hr_CF", "2hr_LTB4", '8hr_CFminus', '8hr_CFplus'))
# drop levels that with no sample
dds$condition <- droplevels(dds$condition)

dds <- DESeq(dds)

resultsNames(dds)

res_8hrp_2hrCF <- results(dds, contrast = c('condition', '8hr_CFplus', "2hr_CF"))
res_8hrp_8hrm <- results(dds, contrast = c('condition', '8hr_CFplus','8hr_CFminus'))

res_list_all_trb <- list('res_8hrp_2hrCF' = res_8hrp_2hrCF  , 'res_8hrp_8hrm' = res_8hrp_8hrm)
#---------------------------------------------------------------------
# just create normalized count
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

## subset gene
gene_name <- read_table(file = 'data/geneset.txt', skip = 2, col_names = 'gene name',
           col_types = cols(
             `gene name` = col_character()
           ))

pos <- rownames(normalized_counts) %in% pull(gene_name)
normalized_counts_sub <- normalized_counts[pos,] 
write.xlsx2(x = normalized_counts_sub, file = 'output/normalized_counts_sub.xlsx')




