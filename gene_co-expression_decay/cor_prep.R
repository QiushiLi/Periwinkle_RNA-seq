##data prep for the pipeline examining relationships between gene expression profiles across different tissues, with Spearman correlation. ##

#####START HERE#####

setwd("/Users/Qiushi/Dropbox/Mac/Documents/croseus/DE/TPM")

library("dplyr")
library("data.table")
library("spaa")

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
