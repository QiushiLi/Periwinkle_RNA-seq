library("dplyr")
combine_TPM <- read.table("/Users/Qiushi/Dropbox/Mac/Documents/croseus/DE/TPM/0906/combine_TPM.txt", sep = "\t")
colnames(combine_TPM) = c('GeneID', 'Root-1', 'Root-2', 'Root-3',	'Leaf-1',	'Leaf-2',	'Leaf-3',	'Epidermis-1', 'Epidermis-2', 'Epidermis-3')
list <- read.table("/Users/Qiushi/Dropbox/Mac/Documents/croseus/DE/TPM/0906/T16H1_MATE-1_CDS.list")
colnames(list) = 'GeneID'
ordered_TPM <- merge(x = list, y = combine_TPM, by = 'GeneID', all.x = T, sort = F)
###add this line to remove duplicated records 
ordered_TPM <- ordered_TPM %>% distinct(GeneID, .keep_all = T)
write.table(ordered_TPM, '/Users/Qiushi/Dropbox/Mac/Documents/croseus/DE/TPM/0906/T16H1_MATE-1_CDS_TPM.txt', sep = '\t', quote = F, row.names = F)
