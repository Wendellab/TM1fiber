library(data.table)
library(tidyverse)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE); # do not omit




groups <- fread("ME.Louvain.Infomap.groups.tsv", header=T, sep="\t")


subgroup <- c("1--10--42",
              "12--4--35",
              "1--4--25",
              "1--4--4",
              "14--4--35",
              "14--9--35",
              "1--6--22",
              "1--6--6",
              "1--9--19",
              "2--1--797",
              "2--2--1583",
              "2--4--712",
              "2--9--1536",
              "2--9--691",
              "2--9--774",
              "3--1--43",
              "4--2--12")

CesGroups <- groups %>% dplyr::filter(group %in% subgroup)


##### ##### ##### ##### ##### ##### ##### ##### 
##### let's deal with ME1 groups first    #####
##### ##### ##### ##### ##### ##### ##### ##### 


ME1 <- fread('cytoscape/CytoscapeInput-edges-ME1',sep="\t",header=T)
CesA <- fread('PCWSCW.tsv')

cesa <- CesA$homoeolog

 
dimensions = function(df) {
  n = c(cesa,df$fromNode,df$toNode) %>% sort %>% unique %>% length
  e = nrow(df)
  paste0(n," nodes, ", e, " edges")
}

write.dimensions = function(df) {
  n = c(df$fromNode,df$toNode) %>% sort %>% unique %>% length
  e = nrow(df)
  fwrite(df, file=paste0("out-12.MEgroup.",n,"nodes.", e, "edges.tsv"), sep="\t")
  paste0(n," nodes, ", e, " edges")
  
}

ME1a <- left_join(ME1,CesA,by=c("fromNode" = "homoeolog")) %>%
  select(-toAltName,-fromAltName) %>%
  rename(fromAltName = CesA, fromStage = stage) %>%
  left_join(.,CesA,by=c("toNode" = "homoeolog")) %>%
  rename(toAltName = CesA, toStage = stage)

onesixsix <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--6--6")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--6--6")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

oneten42 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--10--42")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--10--42")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twelvefour35 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="12--4--35")$homoeolog ) | toNode %in% (filter(CesGroups, group=="12--4--35")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

onefour25 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--4--25")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--4--25")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

onefour4 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--4--4")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--4--4")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

fourteenfour35 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="14--4--35")$homoeolog ) | toNode %in% (filter(CesGroups, group=="14--4--35")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

fourteennine35 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="14--9--35")$homoeolog ) | toNode %in% (filter(CesGroups, group=="14--9--35")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

onesix22 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--6--22")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

onenine19 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="1--9--19")$homoeolog ) | toNode %in% (filter(CesGroups, group=="1--9--19")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twoone797 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--1--797")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--1--797")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twotwo1583 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--2--1583")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--2--1583")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twofour712 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--4--712")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--4--712")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twonine1536 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--9--1536")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--9--1536")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twonine691 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--9--691")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--9--691")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

twonine774 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="2--9--774")$homoeolog ) | toNode %in% (filter(CesGroups, group=="2--9--774")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

threeone43 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="3--1--43")$homoeolog ) | toNode %in% (filter(CesGroups, group=="3--1--43")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))

fourtwo12 <- ME1a %>%
  filter(fromNode %in% (filter(CesGroups, group=="4--2--12")$homoeolog ) | toNode %in% (filter(CesGroups, group=="4--2--12")$homoeolog )) %>%
  filter(weight > quantile(.$weight, probs=0.995))



write.dimensions(onesixsix)
write.dimensions(oneten42)
write.dimensions(twelvefour35)
write.dimensions(onefour25)
write.dimensions(onefour4)
write.dimensions(fourteenfour35)
write.dimensions(fourteennine35)
write.dimensions(onesix22)
write.dimensions(onenine19)
write.dimensions(twoone797)
write.dimensions(twotwo1583)
write.dimensions(twofour712)
write.dimensions(twonine1536)
write.dimensions(twonine691)
write.dimensions(twonine774)
write.dimensions(threeone43)
write.dimensions(fourtwo12)







################ what about using the seidr network

dimensionsdr = function(df) {
  n = c(df$source,df$target) %>% sort %>% unique %>% length
  e = nrow(df)
  paste0(n," nodes, ", e, " edges")
}

write.sdrdimensions = function(df,name) {
  n = c(df$source,df$target) %>% sort %>% unique %>% length
  e = nrow(df)
  nnodes = filter(nodes, nodes$`shared name` %in% (c(df$source,df$target) %>% sort %>% unique))
  fwrite(df, file=paste0("out-12.sdrDirectedgroup.",name,".",n,"nodes.", e, "edges.tsv"), sep="\t")
  paste0(n," nodes, ", e, " edges")
  fwrite(nodes, file=paste0("out-12.sdrDirectedgroup.",name,".",n,"nodes.", e, "edges.NodesOnly.tsv"), sep="\t")
}



sdr <- fread("seidr.top.10pct.tsv")
nodes <- fread("NodeNames.tsv")

#separate scores from ranks
sdr <- sdr %>%
  dplyr::select(source,target,directed,IRP_score) %>%
  left_join(.,CesA,by=c("source" = "homoeolog")) %>% 
  rename(fromAltName = CesA, fromStage = stage) %>%
  left_join(.,CesA,by=c("target" = "homoeolog")) %>%
  rename(toAltName = CesA, toStage = stage)


Sonesix22 <- sdr %>%
  filter(source %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | target %in% (filter(CesGroups, group=="1--6--22")$homoeolog )) %>%
  filter(directed=="Directed") %>%
  filter(IRP_score > quantile(.$IRP_score, probs=0.9))

write.sdrdimensions(Sonesix22, "Sonesix22")

sdrDirected <- sdr %>%
  filter(directed == "Directed")

write.sdrdimensions(sdrDirected, "Sdirectedall")

Sonesix6or22 <- sdr %>%
  filter(source %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | target %in% (filter(CesGroups, group=="1--6--22")$homoeolog )  | 
           target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )  | target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )) %>%
  filter(directed=="Directed") %>%
  filter(IRP_score > quantile(.$IRP_score, probs=0.9))

write.sdrdimensions(Sonesix6or22, "Sonesix6or22")

Sonesix6or22r <- sdr %>%
  filter(source %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | target %in% (filter(CesGroups, group=="1--6--22")$homoeolog )  | 
           target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )  | target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )) %>%
  filter(directed=="Directed") %>%
  filter(IRP_score > quantile(.$IRP_score, probs=0.9))

dimensionsdr(Sonesix6or22r)

write.sdrdimensions(Sonesix6or22, "Sonesix6or22all")



Sonesix6or22u <- sdr %>%
  filter(source %in% (filter(CesGroups, group=="1--6--22")$homoeolog ) | target %in% (filter(CesGroups, group=="1--6--22")$homoeolog )  | 
           target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )  | target %in% (filter(CesGroups, group=="1--6--6")$homoeolog )) %>%
  filter(IRP_score > quantile(.$IRP_score, probs=0.9))

write.sdrdimensions(Sonesix6or22u, "Sonesix6or22undirectedtoo")





# Assuming your dataframe is named 'your_dataframe'
# and it has 75000 rows

# Create a vector of indices to split the dataframe
indices <- seq(1, nrow(nodes), by = 3500)

# Split the dataframe into a list of smaller dataframes
nodechunks <- split(nodes, findInterval(1:nrow(nodes), indices, rightmost.closed = TRUE))

# Save each chunk as a separate file
for (i in seq_along(nodechunks)) {
  filename <- paste0("nodechunk_", i, ".tsv")
  write.table(nodechunks[[i]], file = filename, row.names = FALSE, sep="\t", quote=F)
}


