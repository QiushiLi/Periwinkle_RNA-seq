##data prep for the pipeline examining relationships between gene expression profiles across different tissues, examine the spatial and correlation decay of expression signals##

#####START HERE#####

setwd("/Users/Qiushi/Dropbox/Mac/Documents/croseus/plot/genome_expression/gene_co-expression_decay_croseus")

library("dplyr")
library("data.table")
library("spaa")

################# expression correlation ##############

exp <- read.table("combine_TPM.txt", header=F, sep="\t", row.names=1)
sample_names <- c("Root-1", "Root-2", "Root-3", "Leaf-1", "Leaf-2", "Leaf-3", "Epidermis-1", "Epidermis-2", "Epidermis-3")
setnames(exp, sample_names)

exp <- mutate_all(exp, function(x) as.numeric(as.character(x)))
str(exp)

exp_filt = exp[rowSums(exp >= 0.1) >= 3, ]
dim(exp_filt)

exp_filt_log <- log(exp_filt+1)

cor <- cor(t(exp_filt_log), method="spearman")
corlist <- dist2list(as.dist(cor))

dim(corlist)
head(corlist)

################# gene distance ##############

genes <- fread("cro_v2_cds_to_peri_mRNA.gff3", header = F, sep = "\t")
new_names <- c("chr","source","type","start","end","score","strand","phase","attributes")
setnames(genes, new_names)
genes <- genes[type == "mRNA"]
genes$att <- strsplit(genes$attributes, ";")
genes$geneID <- sapply(genes$att, "[", 3)
genes$geneID <- gsub('Parent=', '', genes$geneID)

genes$gene_TSS = ifelse(genes$strand == "+", genes$start, genes$end)

genes1 <- data.frame(genes$geneID, genes$chr, genes$gene_TSS)
colnames(genes1) <- c("geneID", "chr", "TSS")

splitgenes1 <- split(genes1, genes1$chr)
#test_df = splitgenes1[c('croseus_cScafs_19', 'croseus_cScafs_23')]
distlist_all <- data.frame(col=as.Date(character()),
                           row=character(), 
                           value=character(), 
                           stringsAsFactors=FALSE) 
for (split in splitgenes1){
  mdist <- dist(split$TSS)
  mdist <- as.matrix(mdist)
  rownames(mdist) <- split$geneID
  colnames(mdist) <- split$geneID
  distlist <- dist2list(as.dist(mdist))
  # distmat <- list2dist(mdist)
  distlist_all = rbind(distlist_all, distlist)
}
dim(distlist_all)
head(distlist_all)

mergelist <- merge(distlist_all,corlist, by=c("col", "row"))

write.table(mergelist, "mergelist.txt", sep="\t", quote=FALSE)
save.image(file = "decay.all.test.RData")