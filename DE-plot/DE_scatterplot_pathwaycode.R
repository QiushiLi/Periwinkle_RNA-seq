##used to plot the gene differential expression analysis result for Epidermis, Leaf, and Root tissues of C. roseus, divided into four quadrants, targeted genes are color coded by different pathways##

#####START HERE#####

setwd("/Users/Qiushi/Documents/croseus/DE/job2-franke")

library("ggplot2")
library("DESeq2")
library("dplyr")

###read count_matrix, remember to set sep = ',', when it's a csv file. count_data <- as.matrix(input_data)
input_data <- read.table("gene_count_matrix.csv", header = TRUE, row.names = 1, sep = ',')
count_data <- as.matrix(input_data)

###read samplegroup information
sampleGroup <- read.table("samplegroup.csv", header = TRUE, sep = ",")

### add level information to $Group
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("epidermis","root","leaf"))

### construct dds with DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sampleGroup, design = ~ Group)

### only keep those with ave count > 1 in each replica
keep <- rowSums(counts(dds)) > 3
dds <- dds[keep, ]
dds <- DESeq(dds)

###get comparison results
leaf_vs_root <- results(dds, contrast=c("Group","leaf","root"))
leaf_vs_epidermis <- results(dds, contrast=c("Group","leaf","epidermis"))
epidermis_vs_root <- results(dds, contrast=c("Group","epidermis","root"))

### add row names for the DE result
leaf_vs_root$genes <- rownames(leaf_vs_root)
leaf_vs_epidermis$genes  <- rownames(leaf_vs_epidermis)
epidermis_vs_root$genes  <- rownames(epidermis_vs_root)

### write a separate log2FoldChange dataframe for each comparison
# lgfc_leaf2root <- data.frame(genes=row.names(dds), leaf_vs_root_log2fc=leaf_vs_root$log2FoldChange)
# lgfc_leaf2epi <- data.frame(genes=row.names(dds), leaf_vs_epidermis_log2fc=leaf_vs_epidermis$log2FoldChange)
# lgfc_epi2root <- data.frame(genes=row.names(dds), epidermis_vs_root_log2fc=epidermis_vs_root$log2FoldChange)

###combine lgfc_root, lgfc_epi and epidermis_vs_root with $gene as pivot
# comb_res <- lgfc_leaf2root
# comb_res$leaf_vs_epidermis_log2fc <- lgfc_leaf2epi$leaf_vs_epidermis_log2fc[match(comb_res$genes,lgfc_leaf2epi$genes)]
# comb_res$epidermis_vs_root_log2fc <- lgfc_epi2root$epidermis_vs_root_log2fc[match(comb_res$genes,lgfc_epi2root$genes)]
# rownames(comb_res) = comb_res$genes

###first attempt to plot with separated lists
# list1 <- rownames(subset(leaf_vs_root, padj<0.05 & abs(log2FoldChange)>=1))
# list4 <- rownames(subset(leaf_vs_epidermis, padj<0.05 & abs(log2FoldChange)>=1))

# blabla_plot = ggplot(data=subset(comb_res, genes %in% intersect(list1, list4)), 
#               aes(x=leaf_vs_root_log2fc, y=leaf_vs_epidermis_log2fc)) + 
#               geom_point()
# blabla_plot

### two merge operations to combine 'genes' 'log2FoldChange' 'padj' columns for all comparisons.
combine_results = merge(x=as.data.frame(leaf_vs_root) %>% select(c('genes', 'log2FoldChange', 'padj')),
             y=as.data.frame(leaf_vs_epidermis) %>% select(c('genes', 'log2FoldChange', 'padj')),
             by = c('genes', 'genes'), all = T)

combine_results = merge(x=combine_results,
             y=as.data.frame(epidermis_vs_root) %>% select(c('genes', 'log2FoldChange', 'padj')),
             by = c('genes', 'genes'), all = T)

###add column names for combined results, names() equvalent to colnames()
# names(combine_results) = c('genes', 'lgfc_leaf2root', 'padj_leaf2root', 
                'lgfc_leaf2epi', 'padj_leaf2epi', 'lgfc_epi2root', 'padj_epi2root')
colnames(combine_results) = c('genes', 'lgfc_leaf2root', 'padj_leaf2root', 
                'lgfc_leaf2epi', 'padj_leaf2epi', 'lgfc_epi2root', 'padj_epi2root')

leaf2root_leaf2epi_1_plot = ggplot(data=subset(combine_results, padj_leaf2root<0.05 & abs(lgfc_leaf2root)>=1 &
                                            padj_leaf2epi<0.05 & abs(lgfc_leaf2epi)>=1), 
                              aes(x=lgfc_leaf2root, y=lgfc_leaf2epi)) + geom_point() + labs(title = 'leaf2root_leaf2epi_1')
leaf2root_leaf2epi_1_plot

leaf2root_leaf2epi_2_plot = ggplot(data=subset(combine_results, padj_leaf2root<0.05 & abs(lgfc_leaf2root)>=2 &
                                                 padj_leaf2epi<0.05 & abs(lgfc_leaf2epi)>=2), 
                                   aes(x=lgfc_leaf2root, y=lgfc_leaf2epi)) + geom_point() + labs(title = 'leaf2root_leaf2epi_2')
leaf2root_leaf2epi_2_plot

leaf2root_leaf2epi_epi2root_1 = ggplot(data=subset(combine_results, padj_leaf2root<0.05 & abs(lgfc_leaf2root)>=1 &
                                                     padj_leaf2epi<0.05 & abs(lgfc_leaf2epi)>=1 & 
                                                     padj_epi2root<0.05 & abs(lgfc_epi2root)>=1), 
                                       aes(x=lgfc_leaf2root, y=lgfc_leaf2epi)) + geom_point() + labs(title = 'leaf2root_leaf2epi_epi2root_1')
leaf2root_leaf2epi_epi2root_1

# pathway <- read.table("full1.list", header = F, col.names = 'genes')
otherpw <- read.table("other_pathway.list", header = F, col.names = 'genes')
MEP <- read.table("MEP_gene.list", header = F, col.names = 'genes')
iridoid <- read.table("iridoid_gene.list", header = F, col.names = 'genes')
alkaloid <- read.table("alkaloid_gene.list", header = F, col.names = 'genes')
late_alkaloid <- read.table("late_alkaloid_gene.list", header = F, col.names = 'genes')
ps <- read.table("photosystem_gene.list", header = F, col.names = 'genes')

otherpw <- data.frame(genes=otherpw, cat = 'other_pathway')
MEP <- data.frame(genes=MEP, cat = 'MEP')
iridoid <- data.frame(genes=iridoid, cat = 'iridoid')
alkaloid <- data.frame(genes=alkaloid, cat = 'alkaloid')
late_alkaloid <- data.frame(genes=late_alkaloid, cat = 'late_alkaloid')
ps <- data.frame(genes=ps, cat = 'photosystem')

cat_names <- rbind(otherpw, MEP, iridoid, alkaloid, late_alkaloid, ps)

test <- combine_results

test$cat <- cat_names$cat[match(test$genes, cat_names$genes)]

test_nobackground <- test[!(is.na(test$cat)), ]
# test$cat[is.na(test$cat)] = 'NA'

testcolor = ggplot(data=subset(test_nobackground, padj_leaf2root<0.05 & abs(lgfc_leaf2root)>=0 &
                                 padj_leaf2epi<0.05 & abs(lgfc_leaf2epi)>=0), 
                   aes(x=lgfc_leaf2root, y=lgfc_leaf2epi, color=cat)) + 
  geom_point(shape = 20, alpha = 0.8, size = 2) + labs(title = 'leaf2root_leaf2epi_epi2root_1') + 
  scale_color_manual(name = "cat",
                   values = c("MEP" = "green",
                              "iridoid" = "blue",
                              "alkaloid" = "red",
                              "late_alkaloid" = "brown",
                              "other_pathway" = "grey80",
                              "photosystem" = "darkgreen" ),
                   labels = c("MEP", "iridoid", "alkaloid", "late_alkaloid", "other_pathway", "photosystem")) +
  theme_bw()

testcolor