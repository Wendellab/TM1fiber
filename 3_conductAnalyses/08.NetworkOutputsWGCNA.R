library(RColorBrewer)
library(tidyverse)
library(topGO)
library(Rgraphviz)


options(scipen = 999)
set.seed(8675309)
options(stringsAsFactors = FALSE); # do not omit

# setwd("../filteredNetwork25DPA")
# setwd("W:/corrinne/fiberDan/RNAseq/pseudoAD1/DEanalysis/networks/filteredNetwork25DPA")

############################
## 5. network outputs     ##
############################

load('D5annotation.Rdata')
# "annotation221" "annot"

load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"  "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"

load(file = "04.R-01-dataInput.filtered.RData")
load(file = "06.R-03-buildNetwork.RData")
load(file = "07.R-04-networkTopology.RData")



### Annotate modules with topGO ###

universe <- colnames(datExpr)

# for each module containing genes of interest
GOresults <- data.frame()

# modify geneID2GO
a <- geneID2GO
d <- geneID2GO
names(a) <- paste0(names(a),".1.A")
names(d) <- paste0(names(d),".1.D")
geneID2GOp <- c(a,d)

for(module in 1:(Nmodules-1))
{
    genes <- universe[net$colors==module]
    geneList <- factor(as.integer(universe %in% genes))
    names(geneList) <- universe
    
    pdf(file=paste("topGO/ME",module,".topGO.pdf", sep=""))
    
    # topGO analysis
    remove(enrich)
    
    for(on in c("MF","BP","CC"))
    {
        print(on)
        # Make topGO object
        GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GOp)

        # fisher test
        result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        results.table <- GenTable(GOdata, result, topNodes = length(result@score))

        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
        results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")

        # label ontology type
        results.table$ontology<-on

        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
        keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,]
        if(exists("enrich")) enrich<- rbind(enrich, keep)
        if(!exists("enrich")) enrich<- keep
        
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
    dev.off()

    if(dim(enrich)[1]>0)
    {
        enrichME <- enrich
        enrichME$ME = module
        GOresults <- rbind(GOresults, enrichME)   }
}

GOresults <- GOresults %>% na.omit()

write.table(GOresults, file="out-08.s6.topGO.enrichment.filtered.tsv", sep="\t", row.names=FALSE, quote=F)

parseGO <- GOresults %>% 
  mutate(ME=paste0("ME",ME), GO.ID=paste(GO.ID,qval.bh)) %>% 
  dplyr::select(GO.ID,ME) %>% 
  rowid_to_column(var="rowid") %>% 
  pivot_wider(names_from=ME,values_from=GO.ID) %>% 
  dplyr::select(!rowid) %>% 
  mutate_all(~sort(., na.last = TRUE)) %>% 
  mutate_all(~replace(., is.na(.), ""))

write.table(parseGO, file="out-08.s6.GOterms.tsv", sep="\t", row.names=FALSE, quote=F)

source("Revigo_MEBuildSource.R")
source("function_plotRevigo.R")

for (object in (grep("_scat$", ls(), value = TRUE))) {
  x <- get(object)
  y <- object %>% gsub("_scat","",.) %>% gsub("[BM][PF]_ME","",.) %>% as.numeric()
  z <- object %>% str_sub(start = 1, end = 2)
  RevigoScatter(x,y,z)
}

for (object in (grep("_tree$", ls(), value = TRUE))) {
  x <- get(object)
  y <- object %>% gsub("_tree","",.) %>% gsub("[BM][PF]_ME","",.) %>% as.numeric()
  z <- object %>% str_sub(start = 1, end = 2)
  RevigoTreeMap(x,y,z)
}

#############################
## 7. export to cytoscape  ##
#############################

dir.create("cytoscape")


TOM = TOMsimilarityFromExpr(datExpr, power = 10)
allGenes <- names(datExpr)

# make sure these are set
mergedColors <- labels2colors(net$colors)
Nmodules <- dim(net$MEs)[2]

# export module networks to cytoscape
for(me in 1:(Nmodules-1)) {
  which.module=paste0("ME",me)
  module.color=labels2colors(me)
  
  inME <- is.finite(match(mergedColors, module.color))
  MEgenes <- allGenes[inME]
  
  METOM <- TOM[inME,inME]
  dimnames(METOM) <- list(MEgenes, MEgenes)
  
  MEcyt <- exportNetworkToCytoscape(METOM,
                                    edgeFile = paste0("cytoscape/CytoscapeInput-edges-",which.module),
                                    nodeFile = paste0("cytoscape/CytoscapeInput-nodes-",which.module),
                                    weighted = TRUE,
                                    threshold = 0.1,
                                    nodeNames = MEgenes,
                                    nodeAttr = mergedColors[inME])
  
}













