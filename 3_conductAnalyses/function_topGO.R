library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(tidyverse)

do_topGO <- function(genes,name) {
  geneList <- factor(as.integer(universe %in% genes))
  names(geneList) <- universe
  
  enrich <- data.frame(GO.ID=character(),
                Term=character(), Annotated=integer(),
                Significant=integer(), Expected=numeric(),
                result1=character(), qval.bh=numeric(),
                ontology=character(), category=character(), 
                stringsAsFactors=F)
  
  pdf(file=paste("topGO/",name,".topGO.pdf", sep=""))
  
  # topGO analysis
  
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
    keep <- results.table[as.numeric(results.table[,"qval.bh"])<0.05,] %>% add_column(category=name)
    GOresults <- bind_rows(GOresults,keep)
    
    # draw figure for GO terms pval<=0.05 before FDR correction
    if(is.na(sigNo<-length(keep$ontology))){next}
    showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
    mtext(on, line=-1)
  }
  dev.off() 
 	assign("GOresults",GOresults,envir = globalenv())	

}