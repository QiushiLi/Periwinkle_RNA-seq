##used to plot the PCA analysis for the RNA-seq profiles of Epidermis, Leaf, and Root tissues of C. roseus, three replicas for each tissue##

#####START HERE#####

###read count_matrix as a matrix
input_data <- read.csv("gene_count_matrix.csv", header = TRUE, row.names = 1)
count_data <- as.matrix(input_data)

###read meta information, set Batch and Tissue as factor, specify levels for Tissue
peri_meta <- read.csv("metadata.csv", header = TRUE, row.names = 1)
peri_meta$Batch <- factor(peri_meta$Batch)
peri_meta$Tissue <- factor(peri_meta$Tissue, levels = c("epidermis","root","leaf"))

###construct dds for pca
dds_pca <- DESeqDataSetFromMatrix(countData = count_data, colData = peri_meta, design = ~ Tissue)

###only keep those with ave count > 1 in each replica
keep <- rowSums(counts(dds_pca)) > 3
dds_pca <- dds_pca[keep, ]
dds_pca <- DESeq(dds_pca)

###regularized-logarithm transformation, rlog, good for plotting
dds_pca_rld <- rlog(dds_pca, blind = FALSE)

###PCA plot by Batch and Tissue
plot1_batch <- plotPCA(dds_pca_rld, intgroup = c("Batch"))
plot1_batch
plot2_tissue <- plotPCA(dds_pca_rld, intgroup = c("Tissue"))
plot2_tissue
