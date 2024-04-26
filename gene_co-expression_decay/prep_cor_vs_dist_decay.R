##used to prepare data, relationships between genes analyzed based on the distance of their transcription start sites##

#####START HERE#####

#install.packages("data.table")
#install.packages("spaa")
setwd("/Users/Qiushi/Dropbox/Mac/Documents/croseus/plot/genome_expression/gene_co-expression_decay_croseus/")
library("data.table")
library("spaa")
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

#write.table(distlist_all, "distlist.all.txt", sep="\t", quote=FALSE)
#write.table(distmat, distmat.all.txt, sep="\t", quote=FALSE)