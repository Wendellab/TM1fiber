library(tidyverse)
source("function_KallistoMerge.R")

### use the functions in function_KallistoMerge.R to aggregate files into dataframes

countdata <- countmerge("DEanalysis/counts") 
countdata[,-1] <- round(countdata[,-1],0)

tpmdata <- tpmmerge("DEanalysis/counts") 

write.table(countdata, file="count.Purdue.pseudoAD1.tsv", quote=F, sep="\t", row.names=F)
write.table(tpmdata, file="tpm.Purdue.pseudoAD1.tsv", quote=F, sep="\t", row.names=F)


metadata <- data.frame(id=names(countdata[,-1]), data.frame(id=names(countdata[,-1])) %>% 
	separate(id, c("species", "accession", "plant", "DPA", "rep"))) %>% 
	dplyr::mutate(rep=replace_na(rep,"a")) 


write.table(metadata, file="metadata.Purdue.pseudoAD1.tsv", quote=F, sep="\t", row.names=F)


### record mapped read counts ###
read.counts.millions <- left_join(metadata, 
	(as.data.frame(colSums(countdata[,-1])/1000000) %>% rownames_to_column %>% dplyr::rename(million_read = 2)), 
	by=c("id"="rowname")) 

write.table(read.counts.millions, file="metadata.Purdue.pseudoAD1.withReadCounts.tsv", quote=F, sep="\t", row.names=F)
