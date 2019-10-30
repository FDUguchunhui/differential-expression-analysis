options(java.parameters = "-Xmx10048m")

library(DESeq2)
library(tidyverse)
library(xlsx)
library(magrittr)

#---------------------------------------------------------------------
dat <- read_csv(file = 'data/RNASeq_Trblockade.csv',
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

# 
# write.xlsx(x = res_subgroup(res_8hrp_8hrm,reg_dir = 'up'), file = 'output/res_8hrp_8hrm_up.xlsx')
# write.xlsx(x = res_subgroup(res_8hrp_8hrm,reg_dir = 'down'), file = 'output/res_8hrp_8hrm_down.xlsx')
# write.xlsx(x = res_subgroup(res_8hrp_2hrCF,reg_dir = 'up'), file = 'output/res_8hrp_2hrCF_up.xlsx')
# write.xlsx(x = res_subgroup(res_8hrp_2hrCF,reg_dir = 'down'), file = 'output/res_8hrp_2hrCF_down.xlsx')
# 

#just double check res_subgroup
# pos_1 <- (res_8hrp_8hrm$log2FoldChange > 1)
# pos_2 <- (res_8hrp_8hrm$log2FoldChange < -1)
# sig <- (res_8hrp_8hrm$padj < 0.1)
# sig[is.na(sig)] <- F
# 
# res_8hrp_8hrm[pos_1 & sig, ]
# res_8hrp_8hrm[pos_2 & sig, ]
# 
# res_subgroup(res_8hrp_8hrm, reg_dir = 'down')


res_list_all <- list('res_8hrp_2hrCF' = res_8hrp_2hrCF  , 'res_8hrp_8hrm' = res_8hrp_8hrm)
#---------------------------------------------------------------------

write.xlsx(x = res_8hrp_2hrCF, file = 'output/result_amanitin.xlsx', sheet = 'res_8hrp_2hrCF')
write.xlsx(x = res_8hrp_8hrm, file = 'output/result_amanitin.xlsx', sheet = 'res_8hrp_8hrm', append = T)