library(WGCNA)
library(RColorBrewer)
library(data.table)
library(DESeq2)
library(tidyverse)
library(ggplot2)

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


#############################
## 3. Network Construction ##
#############################

date()

load(file = "04.R-01-dataInput.filtered.RData")

powers <- c(10)
nGenes <- 70000 # 20000

# fixing namespace error
cor <- WGCNA::cor

net = blockwiseModules(datExpr, 
             checkMissingData = TRUE,
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = nGenes,  # 5000 for 4G memory, 20000 for 16G, 30000 for 32 G
             corType = "pearson",
             power = 10, networkType = "signed", TOMType = "signed",
             saveTOMs = TRUE,
             saveTOMFileBase = "filtered.rld_power10_TOM",
             deepSplit = 2,  # default, known to reasonable
             minModuleSize = min(30, ncol(datExpr)/2 ), # default 20, use 30 for transcriptome
             pamStage = TRUE, pamRespectsDendro = TRUE, # default, known to reasonable
             mergeCutHeight = 0.25, # Threshold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75          
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)


# restore namespace
cor <- stats::cor

mergedColors <- labels2colors(net$colors)

pdf(file = "out-06.s3.GeneClusteringWithModuleColors.filtered.pdf", width = 18, height = 18);

plotDendroAndColors(net$dendrograms[[1]],
	mergedColors[net$blockGenes[[1]]], 
	"Modulecolors", main="Module Colors", 
	dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)


dev.off()

save(net, mergedColors,file = "06.R-03-buildNetwork.RData")

