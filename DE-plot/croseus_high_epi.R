setwd("/Users/Qiushi/Documents/croseus/DE/job2-franke/0906/")
###install.packages("ggrepel")
library("ggplot2")
library("DESeq2")
library("dplyr")
library("ggrepel")

###read count_matrix, remember to set sep = ',', when it's a csv file. count_data <- as.matrix(input_data)
input_data <- read.table("gene_count_matrix.csv", header = TRUE, row.names = 1, sep = ',')
count_data <- as.matrix(input_data)

###read sample group information
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

### store row names of the DE result as column "gene"
leaf_vs_root$genes <- rownames(leaf_vs_root)
leaf_vs_epidermis$genes  <- rownames(leaf_vs_epidermis)
epidermis_vs_root$genes  <- rownames(epidermis_vs_root)

### define a function "mrg_selt_col" to combine 'genes' 'log2FoldChange' 'padj' columns for all comparisons.
mrg_selt_col <- function(x) {
  new_list <- list()
  for (i in 1:length(x)) {
    new_list[[i]] <- as.data.frame(x[i])[,c('genes', 'log2FoldChange', 'padj')]
  }
  Reduce(function(x, y) merge(x, y,  by = 'genes', all=TRUE), new_list)
}

MAO <- list(leaf_vs_root, leaf_vs_epidermis, epidermis_vs_root)
combine_results = mrg_selt_col(MAO)

test1 <- combine_TPM %>% mutate_at(c('Root-1', 'Root-2', 'Root-3',	'Leaf-1',	'Leaf-2',	'Leaf-3',	'Epidermis-1', 'Epidermis-2', 'Epidermis-3'), 
                                   function(x) as.numeric(as.character(x)))
test1$Root_mean <- rowMeans(subset(test1, select = c('Root-1', 'Root-2', 'Root-3'), na.rm = TRUE))
test1$Leaf_mean <- rowMeans(subset(test1, select = c('Leaf-1', 'Leaf-2', 'Leaf-3'), na.rm = TRUE))
test1$Epidermis_mean <- rowMeans(subset(test1, select = c('Epidermis-1', 'Epidermis-2', 'Epidermis-3'), na.rm = TRUE))

gene_pos <- read.table("cro_v2_cds_to_peri_final_sorted_gene_pos.txt", header = F)
colnames(gene_pos) = c('chr_id', 'start', 'end', 'strand', 'GeneID')
gene_pos$lg_Root_mean <- log(test1$Root_mean[match(gene_pos$GeneID, test1$GeneID)])
                           
high_root <- as.data.frame(subset(combine_results, padj_leaf2root <= 0.05 & lgfc_leaf2root <= -1)[,"genes"])
colnames(high_root) = 'Root'

high_root_pos <- merge(x = high_root, y = gene_pos, by.x = 'Root', by.y = 'GeneID', all.x = T, sort =F)

install.packages("GenomicRanges")

high_root_pos_DR <- makeGRangesFromDataFrame(high_root_pos,
                                 keep.extra.columns=T,
                                 ignore.strand=F,
                                 seqinfo=NULL,
                                 seqnames.field="chr_id",
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)


chr_id <- paste(rep("cScaf"), seq(1,8,1), sep = '_')
chr_start <- rep(1,8)
chr_end <- c('76892422', '76820121', '72728253', '71689972', '70440602','66496020', '63584754', '59163259')

croseus_genome <- toGRanges(data.frame(chr = chr_id), start = chr_start, end = chr_end)


kp <- plotKaryotype(genome = croseus_genome, chromosomes="cScaf_7")
kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpPoints(kp, data = high_root_pos_DR, y = high_root_pos_DR$lg_Root_mean, data.panel = 1)




high_epi <- subset(test1, padj_leaf2root<0.05 & lgfc_leaf2root >= 1 &
                     padj_leaf2epi<0.05 & lgfc_leaf2epi <= -1)
high_inleaf <- subset(test1, padj_leaf2root<0.05 & lgfc_leaf2root >= 1 &
                        padj_leaf2epi<0.05 & lgfc_leaf2epi >= 0)
