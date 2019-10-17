# this script is design to used with knetics.R
# please run kinetics.R before running


##################### plot

#make the barchart of 4 time point up-regulated gene number
library(ggplot2)
library(tidyr)
dat <- NULL
dat$time <- as.factor(c('1h', '2h', '4h', '6h'))
dat$down <- sapply(res_list_down, nrow)
dat$up <- sapply(res_list_up, nrow)
dat <- as.data.frame(dat)

ggplot(data = dat, aes(time)) +
  geom_bar(aes(time , weight = down, fill = 'red'), show.legend = FALSE) +
  geom_bar(aes(time , weight = up), show.legend = wFALSE) +
  labs(title = "Down regulated genes over time", x = "time", y = "DEGs")

ggplot(data = dat %>% gather(Variable, reg, -time),
       aes(x = time, y = reg, fill = Variable)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.15)) +
  geom_line(aes(x = time, y = reg, group = Variable, color = Variable),
            stat="identity", show.legend = F) +
  geom_point(show.legend = F) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(title = "Down regulated genes over time", x = "time", y = "DEGs")








##find the top 100 up_reg intersection genes expression of different time point
genes_downreg_1h <- tibble.1h_down[(tibble.1h_down$genes %in% genes_down_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_downreg_2h <- tibble.2h_down[(tibble.2h_down$genes %in% genes_down_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_downreg_4h <- tibble.4h_down[(tibble.4h_down$genes %in% genes_down_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_downreg_6h <- tibble.6h_down[(tibble.6h_down$genes %in% genes_down_intersection), ] %>%
  arrange(desc(log2FoldChange))
write.xlsx2(genes_downreg_1h, 'F:\\repos\\Kinetics\\output\\genes_downreg_1h.xlsx')
write.xlsx2(genes_downreg_2h, 'F:\\repos\\Kinetics\\output\\genes_downreg_2h.xlsx')
write.xlsx2(genes_downreg_4h, 'F:\\repos\\Kinetics\\output\\genes_downreg_4h.xlsx')
write.xlsx2(genes_downreg_6h, 'F:\\repos\\Kinetics\\output\\genes_downreg_6h.xlsx')



genes_upreg_1h <- tibble.1h_up[(tibble.1h_up$genes %in% genes_up_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_upreg_2h <- tibble.2h_up[(tibble.2h_up$genes %in% genes_up_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_upreg_4h <- tibble.4h_up[(tibble.4h_up$genes %in% genes_up_intersection), ] %>%
  arrange(desc(log2FoldChange))
genes_upreg_6h <- tibble.6h_up[(tibble.6h_up$genes %in% genes_up_intersection), ] %>%
  arrange(desc(log2FoldChange))
write.xlsx2(genes_upreg_1h, 'F:\\repos\\Kinetics\\output\\genes_upreg_1h.xlsx')
write.xlsx2(genes_upreg_2h, 'F:\\repos\\Kinetics\\output\\genes_upreg_2h.xlsx')
write.xlsx2(genes_upreg_4h, 'F:\\repos\\Kinetics\\output\\genes_upreg_4h.xlsx')
write.xlsx2(genes_upreg_6h, file = 'F:\\repos\\Kinetics\\output\\genes_upreg_6h.xlsx')






# correlation analysis
proteomics <- read.xlsx2(file = 'Proteomics.xlsx', sheetName = 'Sheet1', stringsAsFactors = F,
                         colClasses = c('character', 'character', rep("double", 9))
) %>%
  as_tibble()
count_10h <- read.xlsx2(file = '10hr_normed_counts.xlsx', sheetName = 'counts',
                        colClasses = c('character', rep("double", 16))) %>%
  as_tibble()
proteomics$GENE.ID <- sapply(proteomics$GENE.ID, sub, pattern = '_HUMAN', replacement = '')
intersc <- intersect(proteomics$GENE.ID, count_10h$ID)
x <- count_10h[count_10h$ID %in% intersc, ] %>%
  select(ID,CFASN_TM_1:CFASN_TM_3)
y <- proteomics[proteomics$GENE.ID %in% intersc,] %>%
  select(GENE.ID,CFASN_1:CFASN_3)

intersc.up <- as.list(read.xlsx2('overlap transcriptome proteome.xlsx', header = F,
                                 sheetIndex = 1, stringsAsFactors = F))



##get the correaltions of transcriptomics and proteinomics of intersection genes
# Initiate data frame of correlation results
correlations <- data.frame(matrix(NA, nrow=length(x$ID), ncol=6))
# Assign the first column as the protein names and their site (B2/B4) subset
correlations[ ,1] <- x[,1]
# Assign column names
colnames(correlations) <- c("gene", "Sample Size", "Pearson CC", "Pearson P-Value",
                            "Spearman Rho", "Spearman P-Value")

for (i in 1:nrow(correlations)) {
  cor.pearson <- cor.test(as.numeric(x[i,][,2:4]), as.numeric(y[i,][,2:4]))
  cor.spearman <- cor.test(as.numeric(x[i,][,2:4]), as.numeric(y[i,][,2:4]), method = 'spearman')
  #Col 2. Sample size
  correlations[i,2] <- as.integer(cor.pearson$parameter) + 2
  #Col 3. pcc
  correlations[i,3] <- as.numeric(cor.pearson$estimate)
  #Col 4. pcc p-value
  correlations[i,4] <- as.numeric(cor.pearson$p.value)
  #Col 6. spearman's rho
  correlations[i,5] <- as.numeric(cor.spearman$estimate )
  #Col 7. spearman's rho p-value
  correlations[i,6] <- as.numeric(cor.spearman$p.value )
}

plot(correlations$`Pearson CC`, pch = 16, cex = 0.5)





#get the correaltions of all_time intersection genes
x <- count_10h[count_10h$ID %in% intersc.up$X1, ] %>%
  select(ID,CFASN_TM_1:CFASN_TM_3)
y <- proteomics[proteomics$GENE.ID %in% intersc.up$X1,] %>%
  select(GENE.ID,CFASN_1:CFASN_3)
# Initiate data frame of correlation results
intersc_up.correlations <- data.frame(matrix(NA, nrow=length(x$ID), ncol=6))
# Assign the first column as the protein names and their site (B2/B4) subset
intersc_up.correlations[ ,1] <- x[,1]
# Assign column names
colnames(intersc_up.correlations) <- c("gene", "Sample Size", "Pearson CC", "Pearson P-Value",
                                       "Spearman Rho", "Spearman P-Value")

for (i in 1:nrow(intersc_up.correlations)) {
  cor.pearson <- cor.test(as.numeric(x[i,][,2:4]), as.numeric(y[i,][,2:4]))
  cor.spearman <- cor.test(as.numeric(x[i,][,2:4]), as.numeric(y[i,][,2:4]),
                           method = 'spearman')
  #Col 2. Sample size
  intersc_up.correlations[i,2] <- as.integer(cor.pearson$parameter) + 2
  #Col 3. pcc
  intersc_up.correlations[i,3] <- as.numeric(cor.pearson$estimate)
  #Col 4. pcc p-value
  intersc_up.correlations[i,4] <- as.numeric(cor.pearson$p.value)
  #Col 6. spearman's rho
  intersc_up.correlations[i,5] <- as.numeric(cor.spearman$estimate )
  #Col 7. spearman's rho p-value
  intersc_up.correlations[i,6] <- as.numeric(cor.spearman$p.value )
}

barplot(intersc_up.correlations[intersc_up.correlations$`Pearson P-Value` < 0.05, ]$`Pearson CC`)

intersc_up.correlations[intersc_up.correlations$`Pearson P-Value` < 0.05, ] %>%
  arrange(`Pearson P-Value`) %>%
  write.xlsx2(file = 'proteinomics-transcriptomic-up-pearson-significant.xlsx')

#plot the bubble plot
# par(mfrow = c(2,2))
# data.hallmark.enrich <- read.xlsx2(file = 'Enrichment figures R Chunhui.xlsx',
#                           sheetName = 'Hallmarks',
#                           colClasses = c('character','numeric','character', rep('numeric', 4))
#                           )
# data.hallmark.enrich$FDR.q.value = -log10(as.numeric(data.hallmark.enrich$FDR.q.value))
# data.hallmark.enrich$Gene.Set.Name <- factor(data.hallmark.enrich$Gene.Set.Name,
#                                   levels = data.hallmark.enrich$Gene.Set.Name[order(data.hallmark.enrich$X..Genes.in.Overlap..k.)])
# ggplot(data = data.hallmark.enrich) +
#   geom_point(aes(x = data.hallmark.enrich$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
#   coord_flip() +
#   labs(title = 'Hallmark', size = '# of Genes in Overlap (k)', color = 'FDR q value') +
#   ylab('') +
#   xlab('Genes') +
#   scale_colour_gradient(low = "red", high = "green")
#
# data.comb.enrich <- read.xlsx2(file = 'Enrichment figures R Chunhui.xlsx',
#                                    sheetName = 'Combination',
#                                    colClasses = c('character','numeric','character', rep('numeric', 4))
# )
# data.comb.enrich$FDR.q.value = -log10(as.numeric(data.comb.enrich$FDR.q.value))
# data.comb.enrich$Gene.Set.Name <- factor(data.comb.enrich$Gene.Set.Name,
#                                              levels = data.comb.enrich$Gene.Set.Name[order(data.comb.enrich$X..Genes.in.Overlap..k.)])
# ggplot(data = data.comb.enrich) +
#   geom_point(aes(x = data.comb.enrich$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
#   coord_flip() +
#   labs(title = 'combination', size = '# of Genes in Overlap (k)', color = 'FDR q value') +
#   ylab('') +
#   xlab('Genes') +
#   scale_colour_gradient(low = "red", high = "green")
#
# data.tf.enrich <- read.xlsx2(file = 'Enrichment figures R Chunhui.xlsx',
#                                    sheetName = 'TF',
#                                    colClasses = c('character','numeric','character', rep('numeric', 4))
# )
# data.tf.enrich$Gene.Set.Name <- paste('                               ', data.tf.enrich$Gene.Set.Name)
# data.tf.enrich$FDR.q.value = -log10(as.numeric(data.tf.enrich$FDR.q.value))
# data.tf.enrich$Gene.Set.Name <- factor(data.tf.enrich$Gene.Set.Name,
#                                              levels = data.tf.enrich$Gene.Set.Name[order(data.tf.enrich$X..Genes.in.Overlap..k.)])
# ggplot(data = data.tf.enrich) +
#   geom_point(aes(x = data.tf.enrich$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
#   coord_flip() +
#   labs(title = 'TF', size = '# of Genes in Overlap (k)', color = 'FDR q value') +
#   ylab('') +
#   xlab('Genes') +
#   scale_colour_gradient(low = "red", high = "green") +
#   scale_size(breaks = c(50, 100, 150, 200, 250)) +
#   theme(legend.position = "right")
#
# # plot heatmap
# gene_intsc.table <- read.xlsx2(file = 'heatmap CF upregulated.xlsx', sheetIndex = 1,
#                                stringsAsFactors = F, header = F)
#
# lfc_ts_expr.tibble <- tibble(gene_name = rownames(res_1h),
#       lfc.1h = res_1h$log2FoldChange,
#       lfc.2h = res_2h$log2FoldChange,
#       lfc.4h = res_4h$log2FoldChange,
#       lfc.6h =res_6h$log2FoldChange)
#
# gene_intsc_heatmap.tibble <- lfc_ts_expr.tibble %>%
#   dplyr::filter(gene_name %in% gene_intsc.table[,1])
# gene_intsc_heatmap.matrix <- as.matrix(gene_intsc_heatmap.tibble[,2:5])
# rownames(gene_intsc_heatmap.matrix) <- as.vector(as.list((gene_intsc_heatmap.tibble[,1]))[[1]])
#
# pheatmap(gene_intsc_heatmap.matrix, cluster_cols = F)
# heatmap.2(gene_intsc_heatmap.matrix, Colv = FALSE)










#plot the bubble plot for top 100 genes
#
biological_process.table <- read.xlsx2(file = 'Bubbleplots Chunhui Top100 Up.xlsx',
                                       sheetName = 'BIological process',
                                       colClasses = c('character','numeric','character', rep('numeric', 4))
)

biological_process.table$FDR.q.value = -log10(as.numeric(biological_process.table$FDR.q.value))
#biological_process.table$Gene.Set.Name <- factor(biological_process.table$Gene.Set.Name,
#                                            levels = biological_process.table$Gene.Set.Name[order(biological_process.table$X..Genes.in.Overlap..k.)])

# w:530 h:360
ggplot(data = biological_process.table) +
  geom_point(aes(x = biological_process.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'Biological process', size = '# of Genes in Overlap', color = '-log10 value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 1),
    fill = guide_legend(order = 0)
  )


#*******************************************************************************
cell_component.table <- read.xlsx2(file = 'Bubbleplots Chunhui Top100 Up.xlsx',
                                   sheetName = 'Cel componet',
                                   colClasses = c('character','numeric','character', rep('numeric', 4))
)

cell_component.table$FDR.q.value = -log10(as.numeric(cell_component.table$FDR.q.value))
#cell_component.table$Gene.Set.Name <- factor(cell_component.table$Gene.Set.Name,
#                                                 levels = cell_component.table$Gene.Set.Name[order(data.hallmark.enrich$X..Genes.in.Overlap..k.)])


# w:560 h:400
ggplot(data = cell_component.table) +
  geom_point(aes(x = cell_component.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'Cell component', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )



molecular_function.table <- read.xlsx2(file = 'Bubbleplots Chunhui Top100 Up.xlsx',
                                       sheetName = 'Molecular Function ',
                                       colClasses = c('character','numeric','character', rep('numeric', 4))
)

molecular_function.table$FDR.q.value = -log10(as.numeric(molecular_function.table$FDR.q.value))
#molecular_function.table$Gene.Set.Name <- factor(molecular_function.table$Gene.Set.Name,
#                                                 levels = molecular_function.table$Gene.Set.Name[order(data.hallmark.enrich$X..Genes.in.Overlap..k.)])

# w:780 h:430
ggplot(data = molecular_function.table) +
  geom_point(aes(x = molecular_function.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'molecular function', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )






# another bubble plot for pattern
pattern1.table <- read.xlsx2(file = 'bubble plot patterns.xlsx',
                             sheetName = 'pattern1',
                             colClasses = c('character','numeric','character', rep('numeric', 4))
)

pattern1.table$FDR.q.value = -log10(as.numeric(pattern1.table$FDR.q.value))
#pattern1.table$Gene.Set.Name <- factor(pattern1.table$Gene.Set.Name,
#                                                 levels = pattern1.table$Gene.Set.Name[order(data.hallmark.enrich$X..Genes.in.Overlap..k.)])

# w:515, h:440
ggplot(data = pattern1.table) +
  geom_point(aes(x = pattern1.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 1', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )


# Pattern 2
pattern2.table <- read.xlsx2(file = 'bubble plot patterns.xlsx',
                             sheetName = 'pattern 2',
                             colClasses = c('character','numeric','character', rep('numeric', 4))
)

pattern2.table$FDR.q.value = -log10(as.numeric(pattern2.table$FDR.q.value))
#pattern2.table$Gene.Set.Name <- factor(pattern2.table$Gene.Set.Name,
#                                       levels = pattern2.table$Gene.Set.Name[order(data.hallmark.enrich$X..Genes.in.Overlap..k.)])

ggplot(data = pattern2.table) +
  geom_point(aes(x = pattern2.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 2', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )


#**************************bubble plot
p1_combination.table <- read.xlsx2(file = 'enrichment share up-regulated genes\\bubbleplot data.xlsx',
                                   sheetName = 'p1_combination',
                                   colClasses = c('character','numeric','character', rep('numeric', 4))
)

p1_combination.table$FDR.q.value = -log10(as.numeric(p1_combination.table$FDR.q.value))
# p1_combination.table$Gene.Set.Name <- factor(p1_combination.table$Gene.Set.Name,
#                                        levels = p1_combination.table$Gene.Set.Name[order(p1_combination.table$FDR.q.value)])


#560 440
ggplot(data = p1_combination.table) +
  geom_point(aes(x = p1_combination.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 1 Combination', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )

###############pattern 1 TFT
p1_TFT.table <- read.xlsx2(file = 'enrichment share up-regulated genes\\bubbleplot data.xlsx',
                           sheetName = 'p1_TFT',
                           colClasses = c('character','numeric','character', rep('numeric', 4))
)

p1_TFT.table$FDR.q.value = -log10(as.numeric(p1_TFT.table$FDR.q.value))
#p1_TFT.table$Gene.Set.Name <- factor(p1_TFT.table$Gene.Set.Name,
#                                             levels = p1_TFT.table$Gene.Set.Name[order(p1_TFT.table$FDR.q.value)])

#390 450
ggplot(data = p1_TFT.table) +
  geom_point(aes(x = as.factor(p1_TFT.table$Gene.Set.Name), y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 1 TFT', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )

####pattern 2 combination
p2_combination.table <- read.xlsx2(file = 'enrichment share up-regulated genes\\bubbleplot data.xlsx',
                                   sheetName = 'p2_combination',
                                   colClasses = c('character','numeric','character', rep('numeric', 4))
)

p2_combination.table$FDR.q.value = -log10(as.numeric(p2_combination.table$FDR.q.value))
# p2_combination.table$Gene.Set.Name <- factor(p2_combination.table$Gene.Set.Name,
#                                      levels = p2_combination.table$Gene.Set.Name[order(p2_combination.table$FDR.q.value)])


# 670 440
ggplot(data = p2_combination.table) +
  geom_point(aes(x = p2_combination.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 2 Combination', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )

##pattern 2 TFT
p2_TFT.table <- read.xlsx2(file = 'enrichment share up-regulated genes\\bubbleplot data.xlsx',
                           sheetName = 'p2_TFT',
                           colClasses = c('character','numeric','character', rep('numeric', 4))
)

p2_TFT.table$FDR.q.value = -log10(as.numeric(p2_TFT.table$FDR.q.value))
#p2_TFT.table$Gene.Set.Name <- factor(p2_TFT.table$Gene.Set.Name,
#                                     levels = p2_TFT.table$Gene.Set.Name[order(p2_TFT.table$FDR.q.value)])


#430 440
ggplot(data = p2_TFT.table) +
  geom_point(aes(x =  p2_TFT.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'pattern 2 TFT', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )

# 5-points combination
five_points_combination.table <- read.xlsx2(file = 'enrichment intersection 5 time points genes\\bubbleplot five points.xlsx',
                                            sheetName = 'Combination',
                                            colClasses = c('character','numeric','character', rep('numeric', 4))
)

five_points_combination.table$FDR.q.value = -log10(as.numeric(five_points_combination.table$FDR.q.value))
#five_points_combination.table$Gene.Set.Name <- factor(five_points_combination.table$Gene.Set.Name,
#                                     levels = five_points_combination.table$Gene.Set.Name[order(five_points_combination.table$FDR.q.value)])

#830 440
ggplot(data = five_points_combination.table) +
  geom_point(aes(x =  five_points_combination.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'Combination', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )


# 5-points TFT
five_points_TFT.table <- read.xlsx2(file = 'enrichment intersection 5 time points genes\\bubbleplot five points.xlsx',
                                    sheetName = 'TFT',
                                    colClasses = c('character','numeric','character', rep('numeric', 4))
)

five_points_TFT.table$FDR.q.value = -log10(as.numeric(five_points_TFT.table$FDR.q.value))
five_points_TFT.table$Gene.Set.Name <- factor(five_points_TFT.table$Gene.Set.Name,
                                              levels = five_points_TFT.table$Gene.Set.Name[order(five_points_TFT.table$FDR.q.value)])

#420 440
ggplot(data = five_points_TFT.table) +
  geom_point(aes(x =  five_points_TFT.table$Gene.Set.Name, y='', size =  X..Genes.in.Overlap..k., color = FDR.q.value)) +
  coord_flip() +
  labs(title = 'TFT', size = '# of Genes in Overlap', color = '-log10 q value') +
  ylab('') +
  xlab('Genes') +
  scale_colour_gradient(low = "lightsteelblue1", high = "blue4") +
  guides(
    color = guide_colorbar(order = 0),
    fill = guide_legend(order = 1)
  )







# get the 5 time points genes for time series
genes_five_times_points.ts <- table_nodup[table_nodup$Geneid %in% intersc_unique_10_CFASN_up_four_timepoints$gene_name,] %>%
  rowwise %>%
  mutate(Blood = mean(c(HD1_Blood, HD2_Blood, HD6_Blood)),
         '1 hour' = mean(c(HD1_1H, HD2_1H, HD6_1H)),
         '2 hours' =  mean(c(HD1_2H, HD2_2H, HD6_2H)),
         '4 hours' =  mean(c(HD1_4H, HD2_4H, HD6_4H)),
         '6 hours' =  mean(c(HD1_6H, HD2_6H, HD6_6H))) %>%
  select(Geneid, 'Blood','1 hour', '2 hours', '4 hours', '6 hours')

gene_expr_10h.10h_mean <- gene_expr_10h.table[gene_expr_10h.table$ID %in% intersc_unique_10_CFASN_up_four_timepoints$gene_name,] %>%
  rowwise %>%
  mutate('10 hours' = mean(c(CFASN_TM_1, CFASN_TM_2, CFASN_TM_3,CFASN_TM_5, CFASN_TM_6))) %>%
  select(ID, '10 hours')

five_time_points.ts <- dplyr::inner_join(x = genes_five_times_points.ts, y = gene_expr_10h.10h_mean,
                                         by = c('Geneid' = 'ID'))

# export file for time series
write.table(five_time_points.ts, file = 'enrichment intersection 5 time points genes\\five_time_points.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

#plot

autoHeatmap::hmplot(count_matrix = Trblock_counts_normed, enrichment_xlsx = 'data/overlap.xlsx')
