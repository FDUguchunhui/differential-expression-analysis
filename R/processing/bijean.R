#this option assigns more RAM for java to enable high-RAM-comsuming processing
# especially in writing large dataset to .xlxs
options(java.parameters = "-Xmx8048m")


# # install necessary package
library(BiocManager)
# if you have install DESeq2, uncomment the following line
# BiocManager::install("DESeq2")
# BiocManager::install("biomaRt")
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(xlsx)
library(pheatmap)
library(tidyverse)
library(magrittr)
# run next line if you need to get some help from DESeq2
# browseVignettes("DESeq2")

# before move on
# CF_satus
# condition
# subject
# dat_column_name





#-----------------------------------------------------------------------------------
# read in the dataset
# ignore the warning
dat <- read_csv('data/BIJEAN RNA SEQ COUNTS.csv', col_types = cols(
  .default = col_double(),
  Geneid = col_character(),
  X70 = col_logical()
))

# remove the blank column 
dat <- dat %>% select(-X70)


#********************************************************************************************** 
# # translate ensembl_gene_id
# # check available bioMart databases

listMarts()
ensembl=useMart("ensembl")
# # check available datasets
listDatasets(ensembl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# # check available filters
listFilters(mart)
searchFilters(mart, pattern = 'name')
#
gene_name <- getBM(filters= "ensembl_gene_id",
                 attributes= c("ensembl_gene_id", 'external_gene_name'),
                 values=dat$Geneid,mart= mart)

save(gene_name, file = 'temp/gene_name.rda')
#*******************************************************************************************

load(file = 'temp/gene_name.rda')

# table_merge <- sqldf::sqldf('select *
#              from dat
#              left join gene_name on dat.Geneid = gene_name.ensembl_gene_id')

# merge the those two dataset by ensembl_gene_id
dat <- left_join(x = dat, y = gene_name, by = c('Geneid' = 'ensembl_gene_id'))

dat$external_gene_name <-
  ifelse(is.na(dat$external_gene_name),dat$Geneid,
         dat$external_gene_name)

dat %<>% select(external_gene_name, everything(), -Geneid)

# take the mean of observation have same external_gene_id
# this could take about 1 min
# run the code within * if it is your first time to run it
#******************************************************************************
dat <- dat %>% group_by(external_gene_name) %>% summarise_all(.funs = mean)
save(dat, file = 'temp/dat.rda')
#******************************************************************************
load('temp/dat.rda')


#check whether there is duplicate now, if there is no duplicate, a True should be returned
all(!duplicated(dat[,1]))
#check the first 5 rows
head(dat)

# sort column names to make sure samples from the same subject and condition are close to
#each other
column_order <- sort(colnames(dat))
dat %<>% select(external_gene_name,column_order)

# create count matrix
# have no idea why dimnames parameter doesn't work for tibble
count_matrix <-  as.matrix(dat[, -1])
rownames(count_matrix) <-  pull(dat[,1])
mode(count_matrix) <- 'integer'

# extrace condition and subject level information of study design
dat_column_name <- colnames(dat[-1])
subject <- regmatches(dat_column_name, 
                       m=regexpr(pattern = '.+[^_]+(?=_*[1-9]+)', perl = T, text = dat_column_name))

# exact design information of subject CF status, if subject is CF patient then = CF, otherwise HC, 
#which is shortcut of health control
match = regexpr(pattern = '^CF', perl = T, text = subject)
subject_CF_status <- rep('HC', length(subject))
subject_CF_status[match != -1] <- regmatches(subject, m = match)

# exact design information of condition
condition <- rep('NA', length(subject))
condition_factor <- c('RosSep', 'HC_TRANS', 'CF_TRANS', 'M0', 'M1', 'M2', 'M17')
for(i in 1:length(condition_factor)){
  match = regexpr(pattern = condition_factor[i], perl = T, text = subject)
  condition[match != -1] <- regmatches(subject, m = match)
}


# a multi-level model may not feasible for this study or purpose
# I will leave it here for later usage

# 
# coldata <- data.frame(subject_CF_status = subject_CF_status,
#                       condition = condition,
#                       run = dat_column_name,
#                       row.names = colnames(count_matrix))
# #check whether the name of row match col, a True should be returned
# all(rownames(coldata) == colnames(count_matrix))
# 
# dds <- DESeqDataSetFromMatrix(countData = count_matrix,
#                               colData = coldata,
#                               design = ~ subject_CF_status + condition)
# 
# ddsColl <- collapseReplicates(dds, groupby = replicate, dds$run)
# ddsColl$runsCollapsed

# dds <- DESeq(ddsColl)
# resultsNames(dds)
# results(dds, contrast = list(subject_CF_status, condition))


coldata <- data.frame(condition = paste(subject_CF_status, condition, sep = '_'),
                      run = dat_column_name,
                      row.names = colnames(count_matrix))
#check whether the name of row match col, a True should be returned
all(rownames(coldata) == colnames(count_matrix))

# it will give a warning message, just ignore it 
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)

ddsColl <- collapseReplicates(dds, groupby = subject, dds$run)
ddsColl$runsCollapsed

# pre-filtering
# by removing rows in which there are very few reads, we reduce the memory size of the dds data object,
#   and we increase the speed of the transformation and testing functions within DESeq2
keep <- rowSums(counts(ddsColl)) >= 10
ddsColl <- ddsColl[keep,]

# about 1 min
ddsColl <- DESeq(ddsColl)

resultsNames(ddsColl)


# you can get all the paiwise comparison, it takes a long time to run
# designs <- unique(paste(subject_CF_status, condition, sep = '_'))
# results_list <- list()
# for(i in 1:length(designs)){
#   for(j in 1:length(designs)){
#     if(i > j){
#       contrast <- paste('condition', designs[i],'vs', designs[j], sep = '_')
#       results_list[[contrast]] <- results(ddsColl, 
#                                           contrast = c('condition', designs[i], designs[j]))
#     }
#   }
# }


res.CF_RosSep.HC_RosSep <- results(ddsColl, contrast = c('condition', 'CF_RosSep', 'HC_RosSep'))
res.CF_HC_TRANS.HC_HC_TRANS <- results(ddsColl, contrast = c('condition', 'CF_HC_TRANS', 'HC_HC_TRANS'))
res.CF_CF_TRANS.HC_CF_TRANS <- results(ddsColl, contrast = c('condition', 'CF_CF_TRANS', 'HC_CF_TRANS'))

#HC
res.HC_HC_TRANS.HC_RosSep <- results(ddsColl, contrast = c('condition', 'HC_HC_TRANS', 'HC_RosSep'))
res.HC_CF_TRANS.HC_RosSep <- results(ddsColl, contrast = c('condition', 'HC_CF_TRANS', 'HC_RosSep'))

#CF
res.CF_HC_TRANS.CF_RosSep <- results(ddsColl, contrast = c('condition',  'CF_HC_TRANS', 'CF_RosSep'))
res.CF_CF_TRANS.CF_RosSep <- results(ddsColl, contrast = c('condition',  'CF_CF_TRANS', 'CF_RosSep'))

# macrophages
res.CF_M0.HC_M0 <- results(ddsColl, contrast = c('condition', 'CF_M0', 'HC_M0'))
res.CF_M1.HC_M1 <- results(ddsColl, contrast = c('condition', 'CF_M1', 'HC_M1'))
res.CF_M2.HC_M2 <- results(ddsColl, contrast = c('condition', 'CF_M2', 'HC_M2'))
res.CF_M17.HC_M17 <- results(ddsColl, contrast = c('condition', 'CF_M17', 'HC_M17'))

result_list <- c(res.HC_HC_TRANS.HC_RosSep = res.HC_HC_TRANS.HC_RosSep,
                 res.HC_CF_TRANS.HC_RosSe = res.HC_CF_TRANS.HC_RosSep,
                 
                 res.CF_HC_TRANS.CF_RosSep = res.CF_HC_TRANS.CF_RosSep,
                 res.CF_CF_TRANS.CF_RosSep = res.CF_CF_TRANS.CF_RosSep,
                 
                 res.CF_M0.HC_M0 = res.CF_M0.HC_M0,
                 res.CF_M1.HC_M1 = res.CF_M1.HC_M1,
                 res.CF_M2.HC_M2 = res.CF_M2.HC_M2,
                 res.CF_M17.HC_M17 = res.CF_M17.HC_M17
                 )


#----------------------------------------------------------------------------------------
sheetname <- names(result_list)
for(i in 1:length(result_list)){
  write.xlsx(x = result_list[[i]], file = 'output/bijean_result.xlsx', 
             sheetName = sheetname[i],
             append = T)
}


