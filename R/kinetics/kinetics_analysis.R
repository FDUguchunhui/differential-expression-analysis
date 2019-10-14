# this script is design to used with knetics.R
# please run kinetics.R before running
# this script is not designed to run automatically, just for personal, quick and dirty analysis
# My suggestion is that you should check the code, then modify it or write your own

library(tidyverse)
library(utils)

#get all 4 upregluated genes
res_list <-
  list(
    res_1h = res_1h,
    res_2h = res_2h,
    res_4h = res_4h,
    res_6h = res_6h
  )
res_list_up <- lapply(res_list, res_subgroup, reg_dir = 'up')
#get all 4 downregluated genes
res_list_down <- lapply(res_list, res_subgroup, reg_dir = 'down')



# output corresponding excel file
library(xlsx)
# genes that are upregulated
file_up <-
  paste(getwd(), '/', c("1", '2', '4', '6'), '_up.xlsx', sep = "")
for (i in 1:4) {
  write.xlsx2(res_list_up[i], file = file_up[i], row.names = T)
}
# genes that are downregulated
file_down <-
  paste(getwd(), '/', c("1", '2', '4', '6'), '_down.xlsx', sep = "")
for (i in 1:4) {
  write.xlsx2(res_list_down[i], file = file_down[i], row.names = T)
}





# plot the MAplot with the subset result
plotMA(res_list_up[[1]], ylim = c(-3, 3))
plotMA(res_list[[1]], ylim = c(-3, 3))

# # get the length of the dataset of interest
# # following two lines is demo
# length(res[(res$log2FoldChange > 1)  & res_sig_pos ,]$baseMean)
# length(res[(res$log2FoldChange < -1)  & res_sig_pos ,]$baseMean)




# get the name of downregulated genes
# following demo for 1h
write(file = 'F:/repos/Kinetics/output/gene-downregulated_new.txt', 
      x = row.names(res_list_down[[1]]))


#search for a gene whether in down or up regulation
which(row.names(res_list_up[[1]]) == 'ENSG00000102794')
which(row.names(res_list_down[[1]]) == 'ENSG00000102794')
which(row.names(res) == 'ENSG00000102794')



#plot counts of reads for single gene across groups
plotCounts(dds, gene = which.min(res_list[[1]]$padj), intgroup = "condition")



#transform the data into tibble
#do intersection over all time points with up-regulation
tibble.1h_up <-
  as_tibble(as.data.frame(res_list_up[['res_1h']]), rownames = 'genes')
tibble.2h_up <-
  as_tibble(as.data.frame(res_list_up[['res_2h']]), rownames = 'genes')
tibble.4h_up <-
  as_tibble(as.data.frame(res_list_up[['res_4h']]), rownames = 'genes')
tibble.6h_up <-
  as_tibble(as.data.frame(res_list_up[['res_6h']]), rownames = 'genes')
genes_up_intersection <-
  dplyr::intersect(tibble.1h_up$genes, tibble.2h_up$genes) %>%
  dplyr::intersect(tibble.4h_up$genes) %>%
  dplyr::intersect(tibble.6h_up$genes)
length(genes_up_intersection)



genes_up_intersection.table <-
  table_nodup[table_nodup$Geneid %in% genes_up_intersection, ]
write.xlsx2(genes_intersection.table,
            'genes_up_all_time_points_intersection.xlsx')

#do intersection over all time points with down-regulation
tibble.1h_down <-
  as_tibble(as.data.frame(res_list_down[[1]]), rownames = 'genes')
tibble.2h_down <-
  as_tibble(as.data.frame(res_list_down[[2]]), rownames = 'genes')
tibble.4h_down <-
  as_tibble(as.data.frame(res_list_down[[3]]), rownames = 'genes')
tibble.6h_down <-
  as_tibble(as.data.frame(res_list_down[[4]]), rownames = 'genes')
genes_down_intersection <-
  dplyr::intersect(tibble.1h_down$genes, tibble.2h_down$genes) %>%
  dplyr::intersect(tibble.4h_down$genes) %>%
  dplyr::intersect(tibble.6h_down$genes)
length(genes_down_intersection)
write.xlsx2(genes_up_intersection.table,
            'genes_up_all_time_points_intersection.xlsx')


#do intersection over all time points with down-regulation
tibble.1h_down <-
  as_tibble(as.data.frame(res_list_down[[1]]), rownames = 'genes')
tibble.2h_down <-
  as_tibble(as.data.frame(res_list_down[[2]]), rownames = 'genes')
tibble.4h_down <-
  as_tibble(as.data.frame(res_list_down[[3]]), rownames = 'genes')
tibble.6h_down <-
  as_tibble(as.data.frame(res_list_down[[4]]), rownames = 'genes')
genes_down_intersection <-
  dplyr::intersect(tibble.1h_down$genes, tibble.2h_down$genes) %>%
  dplyr::intersect(tibble.4h_down$genes) %>%
  dplyr::intersect(tibble.6h_down$genes)
length(genes_down_intersection)


genes_down_intersection.table <-
  table_nodup[table_nodup$Geneid %in% genes_down_intersection, ]
write.xlsx2(genes_down_intersection.table,
            'genes_down_all_time_points_intersection.xlsx')

#can be used to check whether the biomart give more gene_name than the original dataset
# length(table$Geneid) length(table_nodup$Geneid)
# sum(grepl(table$Geneid, pattern = 'ENSG00*') == TRUE)
# sum(grepl(table_nodup$Geneid, pattern = 'ENSG00*') == TRUE)


overlap_transcriptome_proteomes <-
  read.xlsx2(
    'overlap transcriptome proteome.xlsx',
    sheetName = 'Sheet1',
    header = F,
    as.data.frame = T
  )
length(which(
  genes_up_intersection %in% overlap_transcriptome_proteomes$X1
))


# 2019/10/11 this is a deep analysis of Knetics_All_HD.csv 
# res_list is all the result from DESeq analysis

res_list <-
  list(
    res_1h = res_1h,
    res_2h = res_2h,
    res_4h = res_4h,
    res_6h = res_6h
  )

# get a subset of it

# all up regulated genes, no requirement on alpha level and log fold change
# not used because we can just subset by genes name with all genes dataset 
# res_list_up <-
#   lapply(
#     res_list,
#     res_subgroup,
#     reg_dir = 'up',
#     alpha = 1,
#     reg_LFC = 0
#   )
# 
# res_list_up <-
#   lapply(
#     lapply(res_list_up, as_tibble, rownames = 'gene_id'),
#     dplyr::select,
#     `gene_id`,
#     `log2FoldChange`
#   )

# this is the dataset with all genes up-regulated
# res_up_all <- Reduce(function(x, y)
#   merge(
#     x,
#     y,
#     by = 'gene_id',
#     all = TRUE,
#     no.dups = T
#   ),
#   x = res_list_up)


# the following code block does the same thing but only those significant genes 
# lfc > 1 and alpha < 0.1
# res_list_up_sig <- lapply(res_list, res_subgroup, reg_dir = 'up')
# res_list_up_sig <-
#   lapply(
#     lapply(res_list_up_sig, as_tibble, rownames = 'gene_id'),
#     dplyr::select,
#     `gene_id`,
#     `log2FoldChange`
#   )

# next get those subset
## 1246
intsc_1246_up <-
  Reduce(f = intersect, lapply(res_list_up_sig, `[[`, 'gene_id'))
length(intsc_1246_up)

## 124  not 6
{
  intsc_124_up <- Reduce(f = intersect,
                         lapply(res_list_up_sig[c('res_1h', 'res_2h', 'res_4h')],
                                `[[`, 'gene_id'))
  length(intsc_124_up)
  set_124 <-
    setdiff(intsc_124_up, res_list_up_sig[['res_6h']]$gene_id)
}

## 126 not 4
{
  intsc_126_up <- Reduce(f = intersect,
                         lapply(res_list_up_sig[c('res_1h', 'res_2h', 'res_6h')], 
                                `[[`, 'gene_id'))
  length(intsc_126_up)
  set_126 <-
    setdiff(intsc_126_up, res_list_up_sig[['res_4h']]$gene_id)
}

# 146 not 2
{
  intsc_146_up <-
    Reduce(f = intersect, lapply(res_list_up_sig[c('res_1h', 'res_4h', 'res_6h')], 
                                 `[[`, 'gene_id'))
  length(intsc_146_up)
  set_146 <-
    setdiff(intsc_146_up, res_list_up_sig[['res_2h']]$gene_id)
}


## 246 not 1
{
  intsc_246_up <-
    Reduce(f = intersect, lapply(res_list_up_sig[c('res_2h', 'res_4h', 'res_6h')], 
                                 `[[`, 'gene_id'))
  length(intsc_246_up)
  set_246 <-
    setdiff(intsc_246_up, res_list_up_sig[['res_1h']]$gene_id)
}

# 12
{
  intsc_12_up <-
    Reduce(f = intersect, lapply(res_list_up_sig[c('res_2h', 'res_4h', 'res_6h')], 
                                 `[[`, 'gene_id'))
  length(intsc_12_up)
  set_12 <-
    setdiff(intsc_12_up,
            union(res_list_up_sig[['res_4h']]$gene_id,
                  res_list_up_sig[['res_6h']]$gene_id))
}





combn(c('1', '2', '4', '6'), m = 3)
combn(c('1', '2', '4', '6'), m = 2)



# use for intersection
# res_list <-
#   list(
#     res_1h = res_1h,
#     res_2h = res_2h,
#     res_4h = res_4h,
#     res_6h = res_6h
#   )
# res_list_down <-
#   lapply(
#     res_list,
#     res_subgrodown,
#     reg_dir = 'down',
#     alpha = 1,
#     reg_LFC = 0
#   )
# 
# 
# res_list_down <-
#   lapply(
#     lapply(res_list_down, as_tibble, rownames = 'gene_id'),
#     dplyr::select,
#     `gene_id`,
#     `log2FoldChange`
#   )
# 
# res_down_all <- Reduce(function(x, y)
#   merge(
#     x,
#     y,
#     by = 'gene_id',
#     all = TRUE,
#     no.ddowns = T
#   ),
#   x = res_list_down)
# 
# Reduce(f = intersect, res_list_down)



#
res_list_all <-
  lapply(
    res_list,
    res_subgroup,
    reg_dir = 'all',
    alpha = 1,
    reg_LFC = 0
  )

res_list_all <-
  lapply(
    lapply(res_list_all , as_tibble, rownames = 'gene_id'),
    dplyr::select,
    `gene_id`,
    `log2FoldChange`
  )

res_all <- Reduce(function(x, y)
  merge(
    x,
    y,
    by = 'gene_id',
    all = TRUE,
    no.dups = T
  ),
  x = res_list_all) %>% as_tibble()

colnames(res_all) <- c('gene_id', '1h', '2h', '4h', '6h')



<<<<<<< HEAD
all <- list(res_all[res_all$gene_id %in% intsc_1246_up, ],
            res_all[res_all$gene_id %in% set_124, ],
            res_all[res_all$gene_id %in% set_126, ],
            res_all[res_all$gene_id %in% set_146, ],
            res_all[res_all$gene_id %in% set_246, ])

write.xlsx(all[[1]],  sheetName = 'intsc_1246_up',
           file = 'output/dontknowhowtoname.xlsx')
write.xlsx(all[[2]],
           sheetName = '124',
           file = 'output/dontknowhowtoname.xlsx',
           append = T)
write.xlsx(all[[3]],
           sheetName = '126',
           file = 'output/dontknowhowtoname.xlsx',
           append = T)
write.xlsx(all[[4]],
           sheetName = '146',
           file = 'output/dontknowhowtoname.xlsx',
           append = T)
write.xlsx(all[[5]],
           sheetName = '246',
           file = 'output/dontknowhowtoname.xlsx',
           append = T)
=======
all <- list(res_all[res_all$gene_id %in% intsc_1246_up,],
res_all[res_all$gene_id %in% set_124,],
res_all[res_all$gene_id %in% set_126,],
res_all[res_all$gene_id %in% set_146,],
res_all[res_all$gene_id %in% set_246,]
)

write.xlsx(all[[1]],  sheetName = 'intsc_1246_up', 
                  file = 'output/dontknowhowtoname.xlsx')
write.xlsx(all[[2]],  sheetName = '124', 
           file = 'output/dontknowhowtoname.xlsx', append = T)
write.xlsx(all[[3]],  sheetName = '126', 
           file = 'output/dontknowhowtoname.xlsx', append = T)
write.xlsx(all[[4]],  sheetName = '146', 
           file = 'output/dontknowhowtoname.xlsx', append = T)
write.xlsx(all[[5]],  sheetName = '246', 
           file = 'output/dontknowhowtoname.xlsx', append = T)

#---------------------------------------------------------------------------------------



upreg <- read.xlsx(file = 'data/venn_result Upregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
upreg <- upreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_up <- levels(factor(upreg$Names, levels = unique(upreg$Names)))

up_comb <- list()
for(i in 1:length(combo_up)){
  up_comb[[i]] <- res_all[res_all$gene_id %in% upreg[upreg$Names == combo_up[i],]$elements,]
}


# 
for(i in 1:length(combo_up)){
  write.xlsx(up_comb[[i]],  sheetName = sub('Hours', '', combo_up)[[i]],  # the sheetname is short
             file = 'output/up_comb.xlsx', append = T)
}


# for down reg

downreg <- read.xlsx(file = 'data/Venn_Downregulated.xlsx', sheetIndex = 1, stringsAsFactors = F)
downreg <- downreg %>% as_tibble() %>% fill(`Names`, .direction = c("down"))  
combo_down <- levels(factor(downreg$Names, levels = unique(downreg$Names)))

down_comb <- list()
for(i in 1:length(combo_down)){
  down_comb[[i]] <- res_all[res_all$gene_id %in% downreg[downreg$Names == combo_down[i],]$elements,]
}


# 
for(i in 1:length(combo_down)){
  write.xlsx(down_comb[[i]],  sheetName = gsub('hour[s]*', '', combo_down)[[i]],  # the sheetname is short
             file = 'output/down_comb.xlsx', append = T)
}
>>>>>>> origin/new_method
