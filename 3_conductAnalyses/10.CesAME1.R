library(data.table)
library(tidyverse)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE); # do not omit

ME1 <- fread('cytoscape/CytoscapeInput-edges-ME1',sep="\t",header=T)
CesA <- fread('PCWSCW.tsv')

cesa <- CesA$homoeolog

 
dimensions = function(df) {
  n = c(cesa,df$fromNode,df$toNode) %>% sort %>% unique %>% length
  e = nrow(df)
  paste0(n," nodes, ", e, " edges")
}

ME1a <- left_join(ME1,CesA,by=c("fromNode" = "homoeolog")) %>%
  select(-toAltName,-fromAltName) %>%
  rename(fromAltName = CesA, fromStage = stage) %>%
  left_join(.,CesA,by=c("toNode" = "homoeolog")) %>%
  rename(toAltName = CesA, toStage = stage)
  
# extract the top 5 edges going into CesA
LevelZero <- ME1a %>%
  group_by(toAltName) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))

lvl0 <- c(LevelZero$fromNode, LevelZero$toNode) %>% unique

# extract the top 5 edges going into CesA or a node going directly into it
LevelOne <- ME1a %>%
  filter(toNode %in% lvl0) %>%
  group_by(toNode) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))

lvl1 <- c(LevelOne$fromNode, LevelOne$toNode) %>% unique

# extract the top 5 edges going into CesA or the first two nodes going into it
LevelTwo <- ME1a %>%
  filter(toNode %in% lvl1) %>%
  group_by(toNode) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))



# extract the top 5 edges going out of CesA (because this is not directed, these could be "in" nodes)
LevelZeroB <- ME1a %>%
  group_by(fromAltName) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))

lvl0B <- c(LevelZeroB$fromNode, LevelZeroB$toNode) %>% unique

LevelMinusOne <- ME1a %>%
  filter(toNode %in% lvl0B) %>%
  group_by(toNode) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))

lvlm1 <- c(LevelMinusOne$fromNode, LevelMinusOne$toNode) %>% unique

LevelMinusTwo <- ME1a %>%
  filter(toNode %in% lvlm1) %>%
  group_by(toNode) %>%
  top_n(5, wt = weight) %>%
  ungroup() %>%
  mutate(across(c(fromAltName), ~ifelse(is.na(.), fromNode, .)))
  
  
TwoEachDir <- bind_rows(LevelTwo,LevelMinusTwo) %>% unique

dimensions(LevelZero)      # [1] "95 nodes, 80 edges"
dimensions(LevelOne)       # [1] "188 nodes, 375 edges"
dimensions(LevelTwo)       # [1] "260 nodes, 835 edges"
dimensions(LevelZeroB)     # [1] "87 nodes, 80 edges"
dimensions(LevelMinusOne)  # [1] "204 nodes, 335 edges"
dimensions(LevelMinusTwo)  # [1] "282 nodes, 920 edges"
dimensions(TwoEachDir)     # [1] "319 nodes, 1165 edges"

fwrite(LevelZero,file="ME1.CesA.top5.Lvl0.95n.80e.tsv",sep="\t")
fwrite(LevelOne,file="ME1.CesA.top5.Lvl1.188n.375e.tsv",sep="\t")
fwrite(LevelTwo,file="ME1.CesA.top5.Lvl2.260n.835e.tsv",sep="\t")
fwrite(LevelZeroB,file="ME1.CesA.top5.Lvl0B.87n.80e.tsv",sep="\t")
fwrite(LevelMinusOne,file="ME1.CesA.top5.Lvlm1.204n.335e.tsv",sep="\t")
fwrite(LevelMinusTwo,file="ME1.CesA.top5.Lvlm2.282n.920e.tsv",sep="\t")
fwrite(TwoEachDir,file="ME1.CesA.top5.all.319n.1165e.tsv",sep="\t")


quantile(TwoEachDir$weight, seq(0,1,.1))
#         0%        10%        20%        30%        40%        50%        60%
# 0.02738306 0.22410405 0.26409487 0.30073401 0.32046783 0.33592941 0.34674592
#        70%        80%        90%       100%
# 0.35417107 0.36041925 0.36480682 0.57934998


dimensions(ME1a)
# 
quantile(ME1a$weight, seq(0,1,0.1))
#         0%        10%        20%        30%        40%        50%        60%
# 0.01000000 0.01408644 0.01934317 0.02611281 0.03484655 0.04614056 0.06088163
#        70%        80%        90%       100%
# 0.08064296 0.10881919 0.15408204 0.57934998




####################################
### now redo with only top 3,2,1 ###
####################################

top3 <- TwoEachDir %>%
  group_by(toNode) %>%
  top_n(3, wt = weight) %>%
  ungroup() 
  
top2 <- TwoEachDir %>%
  group_by(toNode) %>%
  top_n(2, wt = weight) %>%
  ungroup() 

top1 <- TwoEachDir %>%
  group_by(toNode) %>%
  top_n(1, wt = weight) %>%
  ungroup() 
  
dimensions(top3) # [1] "278 nodes, 701 edges"
dimensions(top2) # [1] "266 nodes, 466 edges"
dimensions(top1) # [1] "258 nodes, 234 edges"

fwrite(top3,file="ME1.CesA.top3.all.278n.701e.tsv",sep="\t")
fwrite(top2,file="ME1.CesA.top2.all.266n.466e.tsv",sep="\t")
fwrite(top1,file="ME1.CesA.top1.all.258n.234e.tsv",sep="\t")


#######################################
### get a final list of nodes/edges ###
#######################################

MEces2 <- ME1a %>%
	filter(weight > 0.2641) # TwoEachDir 20% quantile
	
MEces3 <- ME1a %>%
	filter(weight > 0.3007)  # TwoEachDir 30% quantile
	
MEces5 <- ME1a %>%
	filter(weight > 0.3359)  # TwoEachDir 50% quantile
	
dimensions(MEces2)  # [1] "9060 nodes, 1412938 edges"
dimensions(MEces3)  # [1] "5834 nodes, 359491 edges"
dimensions(MEces5)  # [1] "2470 nodes, 39775 edges"


fwrite(MEces2,file="ME1.CesA.20thpct.9060n.1412938e.tsv",sep="\t")  
fwrite(MEces3,file="ME1.CesA.30thpct.5834n.359491e.tsv",sep="\t")   
fwrite(MEces5,file="ME1.CesA.50thpct.2470n.39775e.tsv",sep="\t")   
