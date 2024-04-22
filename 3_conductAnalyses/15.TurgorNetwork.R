library(data.table)
library(tidyverse)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE); # do not omit




groups <- fread("ME.Louvain.Infomap.groups.tsv", header=T, sep="\t")


dimensionsdr = function(df) {
  n = c(df$source,df$target) %>% sort %>% unique %>% length
  e = nrow(df)
  paste0(n," nodes, ", e, " edges")
}

write.sdrdimensions = function(df,name) {
  n = c(df$source,df$target) %>% sort %>% unique %>% length
  e = nrow(df)
#  nnodes = filter(nodes, nodes$`shared name` %in% (c(df$source,df$target) %>% sort %>% unique))
  fwrite(df, file=paste0("out-14.sdrDirectedgroup.Turgor.",name,".",n,"nodes.", e, "edges.tsv"), sep="\t")
  paste0(n," nodes, ", e, " edges")
  fwrite(nodes, file=paste0("out-14.sdrDirectedgroup.Turgor.",name,".",n,"nodes.", e, "edges.NodesOnly.tsv"), sep="\t")
}



sdr <- fread("seidr.top.10pct.tsv")
nodes <- fread("NodeNames.tsv")
ME8 <- fread('cytoscape/CytoscapeInput-nodes-ME8',sep="\t",header=T)
ME9 <- fread('cytoscape/CytoscapeInput-nodes-ME9',sep="\t",header=T)
ME89 <- c(ME8$nodeName,ME9$nodeName)

#separate scores from ranks
sdrDir <- sdr %>%
  dplyr::select(source,target,directed,IRP_score) %>%
  filter(directed == "Directed")

sdr89 <- sdrDir %>%
  filter(source %in% ME89 | target %in% ME89)

write.sdrdimensions(sdr89,"ME89")





Sonesix22 <- sdr %>%
  filter(source %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | target %in% (filter(CesGroups, group=="1--6--22")$homoeolog )) %>%
  filter(directed=="Directed") %>%
  filter(IRP_score > quantile(.$IRP_score, probs=0.9))

write.sdrdimensions(Sonesix22, "Sonesix22")

sdrDirected <- sdr %>%
  filter(directed == "Directed")

write.sdrdimensions(sdrDirected, "Sdirectedall")



