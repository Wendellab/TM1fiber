library(WGCNA)
library(RColorBrewer)
library(data.table)
library(DESeq2)
library(tidyverse)

library(doParallel)
registerDoParallel()


options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE); # do not omit
enableWGCNAThreads(nThreads=18); # set WGCNA threads

setwd("../filteredNetwork25DPA")
# setwd("W:/corrinne/fiberDan/RNAseq/pseudoAD1/DEanalysis/networks/filteredNetwork25DPA")


##################
## 1. Check Set ##
##################


### Input read count and coldata ###

count <- fread("count.Purdue.pseudoAD1.filtered.25DPA.tsv", sep="\t")
rownames(count) <- count$target_id

coldata <- fread("metadata.Purdue.pseudoAD1.filtered.25DPA.tsv", sep="\t")
rownames(coldata) <- coldata$id


### count the number of A and D homoeologs ###
idxa <- grep( "[.]A$",rownames(count) )
idxd <- grep( "[.]D$",rownames(count) )


### get the number of total reads and proportion of A vs D ###
coldata$libSize <- colSums(count[,-1])
coldata$AbyD <- colSums(count[idxa,-1])/colSums(count[idxd,-1])

table(coldata$DPA)
coldata[coldata$libSize<5000000,]


### DEseq2 normalization/rlog transformation of count data ###
dds <- DESeqDataSetFromMatrix(countData=count, colData=coldata, design=~1, tidy = TRUE)
rld <- rlog(dds)
count.rld <- as.data.frame(assay(rld))
colnames(count.rld) <- coldata$id

### transform the DEseq normalization to a WGCNA data matrix ###
datExpr <- as.data.frame(t(count.rld))


### Check to see which genes are ok for network construction ###
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK


# Exclude ZERO genes from the calculation due to too many missing samples or zero variance
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
    
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))

    # Remove the offending genes and samples
        datExpr = datExpr[gsg$goodSamples, gsg$goodGenes];
}


gsg = goodSamplesGenes(datExpr, verbose = 3);
str(gsg)
# List of 3
# $ goodGenes  : logi [1:69245] TRUE TRUE TRUE TRUE TRUE TRUE ...
# $ goodSamples: logi [1:55] TRUE TRUE TRUE TRUE TRUE TRUE ...
# $ allOK      : logi TRUE


### Check clustering of samples ###
pdf(file = "out-04.s1.SampleClustering.filtered.pdf", width = 12, height = 12);

sampleTrees = hclust(dist(datExpr), method = "complete")

plot(sampleTrees, 
	main = paste("Sample clustering on all genes"),
	xlab="", sub="", cex = 0.7)
dev.off()


### save output ###
save(datExpr, coldata, count.rld, idxa, idxd, file = "04.R-01-dataInput.filtered.RData")


