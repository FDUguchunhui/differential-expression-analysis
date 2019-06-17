# # install necessary package
library(BiocManager)

# if you have install DESeq2, uncomment the followiing line
# BiocManager::install("DESeq2")
library(DESeq2)

# run next line if you need to get some help from DESeq2
# browseVignettes("DESeq2")

#this option assigns more RAM for java to enable high-RAM-comsuming processing
# especially in writing large dataset to .xlxs
options(java.parameters = "-Xmx2048m")


# read in the dataset
table <- read.csv('kinetics_All_HD.csv', stringsAsFactors = F)


#the original dataset has two gene id, 
#and the second one has duplicates, in order to use the second one, need to deal with it
#!
table_geneid <- table[,-1]
# omit any row if has missing value 
# Caution! this may not the best method
table_geneid <- na.omit(table_geneid)
head(table_geneid)


#remove depublicate by geneid using sqlite
# base rule: average the read from the same gene
# caution! This may not applicable for bio-duplicate
table_nodup <-  sqldf::sqldf('
  select geneid,
         avg(HD1_Blood) as HD1_Blood, avg(HD2_Blood) as HD2_Blood, avg(HD6_Blood) as HD6_Blood,
         avg(HD1_1H) as HD_1H, avg(HD2_1H) as HD2_1H, avg(HD6_1H) as HD6_1H,
         avg(HD1_2H) as HD_2H, avg(HD2_2H) as HD2_2H, avg(HD6_2H) as HD6_2H,
         avg(HD1_4H) as HD_4H, avg(HD2_4H) as HD2_4H, avg(HD6_4H) as HD6_4H,
         avg(HD1_6H) as HD_6H, avg(HD2_6H) as HD2_6H, avg(HD6_6H) as HD6_6H
  from table_geneid
  group by Geneid
  order by Geneid
  '
)

#check how many gene left in the processed dataset
length(table_nodup$Geneid)

#check whether there is duplicate now, if there is no duplicate, a True should be returen
all(!duplicated(table_nodup[,1])) 
#check the first 5 rows
head(table_nodup)






#make the matrix needed for next step

#use the geneid as row.names 
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
# may can be done by lapply() function, not sure how to set the func argument with different arguments
res_1h <- results(dds,contrast=c("condition","1h","blood"))
res_2h <- results(dds,contrast=c("condition","2h","blood"))
res_4h <- results(dds,contrast=c("condition","4h","blood"))
res_6h <- results(dds,contrast=c("condition","6h","blood"))

resultsNames(dds)

# log fold change shrinkage(LFC Shrinkage) may not neccessary here. for using it, uncomment next two lines 
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
  
  return(res[res_LFC_pos,])
  
}




#get all 4 upregluated genes
res_list <- list(res_1h,res_2h,res_4h,res_6h)
res_list_up <- lapply(res_list, res_subgroup, reg_dir ='up')
#get all 4 downregluated genes
res_list_down <- lapply(res_list, res_subgroup, reg_dir ='down')


# output corresponding excel file
library(xlsx)
# genes that are upregulated
file_up <- paste(getwd(),'/', c("1", '2', '4', '6'), '_up.xlsx', sep="")
for(i in 1:4){
  write.xlsx2(res_list_up[i], file = file_up[i], row.names = T)
}
# genes that are downregulated
file_down <- paste(getwd(),'/', c("1", '2', '4', '6'), '_down.xlsx', sep="")
for(i in 1:4){
  write.xlsx2(res_list_down[i], file = file_down[i], row.names = T)
}


# plot the MAplot with the subset result
plotMA(res_list_up[[1]], ylim=c(-3,3))

# get the length of the dataset of interest
# following two lines is demo
length(res[(res$log2FoldChange > 1)  & res_sig_pos ,]$baseMean)
length(res[(res$log2FoldChange < -1)  & res_sig_pos ,]$baseMean)

# get the name of downregulated genes
# following demo for 1h
write(file ='gene-downregulated_new',x = row.names(res_list_down[[1]]))

#search for a gene whether in down or up regulation
which(row.names(res_list_up[[1]]) == 'ENSG00000102794')
which(row.names(res_list_down[[1]]) == 'ENSG00000102794')
which(row.names(res) == 'ENSG00000102794')



#plot counts of reads for single gene across groups
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

