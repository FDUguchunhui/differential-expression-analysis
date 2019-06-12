library(BiocManager)
BiocManager::install("DESeq2")
browseVignettes("DESeq2")

library(DESeq2)
browseVignettes("DESeq2")


table <- read.csv('kinetics_All_HD.csv', stringsAsFactors = F)
#the original dataset has two gene id, 
#and the second one has duplicates, in order to use the second one, need to deal with it
#!
table_geneid <- table[,-1]
table_geneid <- na.omit(table_geneid)
head(table_geneid)


"table_geneid_lite <- table_geneid[1:1000,]
table_nodup <-  table_geneid[0,]
rep <- NULL
for(i in seq(length(table_geneid_lite$Geneid))){
  if(!(table_geneid_lite[i,]$Geneid %in% table_nodup$Geneid)){
    table_nodup[i,] = table_geneid_lite[i,]
  }
  else{
    table_nodup['table_geneid_lite[i]$Geneid',-1] = table_nodup[i,-1] + table_geneid_lite[i,-1]
  }
}"


#remove depublicate by geneid
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

length(table_nodup$Geneid)

#check whether there is duplicate now
all(!duplicated(table_nodup[,1])) 
#check the first 5 rows
head(table_nodup)

#use the geneid1 as row.names (geneid cannot used as row.names because of depulicate)
table_rownames <- data.frame(table_nodup[,-1], row.names=table_nodup[,1])
count_maxtrix <- as.matrix(table_rownames)
mode(count_maxtrix) <- 'integer'
head(count_maxtrix)

#compare blood and 1h data
coldata <- data.frame(condition = c(rep('blood', 3),
                                    rep('1h', 3),
                                    rep('2h', 3),
                                    rep('4h', 3),
                                    rep('6h', 3)),
                      row.names = colnames(count_maxtrix)[1:15])
#check whether the name of row match col
all(rownames(coldata) == colnames(count_maxtrix)[1:15])

#create dds
dds <- DESeqDataSetFromMatrix(countData = count_maxtrix[,1:15],
                              colData = coldata,
                              design = ~ condition)

#pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#set dds condition to all time points
dds$condition <- factor(dds$condition, levels = c("blood","1h",'2h', '4h', '6h'))
dds$condition <- droplevels(dds$condition)

# DEseq analysis
dds <- DESeq(dds)
#for different condition you just change 6h to i.e. 1h, 2h, 4h
res <- results(dds,contrast=c("condition","6h","blood"))
res
resultsNames(dds)


# # log fold change shrinkage
# resLFC <- lfcShrink(dds, coef="condition_1h_vs_blood", type="apeglm")
# resLFC


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




#get those LFC greate than 1 and the adjusted P-value is less than 0.1
res_sig_pos <- (res$padj < 0.1)
res_sig_pos[is.na(res_sig_pos)] <- F
res_LFC_pos <- (res$log2FoldChange > 1) | (res$log2FoldChange < -1)
plotMA(res[ res_LFC_pos & res_sig_pos ,], ylim=c(-3,3))

#upregulated and downregulated with adj P < 0.1
length(res[(res$log2FoldChange > 1)  & res_sig_pos ,]$baseMean)
length(res[(res$log2FoldChange < -1)  & res_sig_pos ,]$baseMean)

#get the name of downregulated genes
downreg_gene <- rownames(res[(res$log2FoldChange < -1)  & res_sig_pos ,])
write(file ='gene-downregulated_new',x = downreg_gene )








#plot counts of reads for single gene across groups
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

