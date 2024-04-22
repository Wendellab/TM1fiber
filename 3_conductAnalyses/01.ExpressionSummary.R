library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)

options(scipen = 999)
set.seed(8675309)
setDTthreads(12) 
options(datatable.fread.datatable=FALSE)
max.overlaps=Inf




tpm <- fread("tpm.Purdue.pseudoAD1.tsv", header=T)
names(tpm) <- gsub("_tpm","",names(tpm))

meta <- fread("metadata.Purdue.pseudoAD1.tsv", header=T)

tbl <- meta %>% 
	select(c(id,DPA)) %>%
	arrange(DPA,id)

withZeroes <- tpm %>% 
  select(-target_id) %>% 
  summarise_all(funs(min = min(.), 
                               q25 = quantile(., 0.25), 
                               median = median(.), 
                               medianGT0 = median(.[.>0]), 
                               medianq5.95 = median(.[.> quantile(., 0.05) & .< quantile(., 0.95)]), 
                               q75 = quantile(., 0.75), 
                               max = max(.), 
                               mean = mean(.), 
                               meanGT0 = mean(.[.>0]),
                               meanq5.95 = mean(.[.> quantile(., 0.05) & .< quantile(., 0.95)]), 
                               sd = sd(.), 
                               q5 = quantile(., 0.05), 
                               q95 = quantile(., 0.95)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=182),
	rep(LETTERS[1:2],times=13),
	rep(LETTERS[1:3],times=39),
	rep(LETTERS[1:2],times=13),
	rep(LETTERS[1:1],times=13) ))

summaryPlot <- ggplot() + 
  geom_line(data=subset(withZeroes,statistic %in% "mean"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "median"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dashed" ) + 
  geom_line(data=subset(withZeroes,statistic %in% "q5"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "q25"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "q75"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dotted" ) + 
  geom_line(data=subset(withZeroes,statistic %in% "q95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "meanGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "medianGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "meanq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroes,statistic %in% "medianq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  scale_color_brewer(palette = "Dark2")

ggsave(summaryPlot, file="out-01.TPM_summaryPlot.jpg", width=6, height=3, units="in")

summaryTable <- withZeroes %>% 
  mutate(sample=paste0(DPA,"_",rep)) %>%
  select(sample,statistic,V1) %>%
  spread(sample,V1)

write.table(file="out-01.summaryTable.tsv", summaryTable, quote=F, row.names=F, sep="\t")



tpm %>% 
  select(-target_id) %>% 
  summarise_all(funs(q5GT0 = quantile(.[.>0],0.05)), 
	            	 na.rm = TRUE) %>% 
     t() %>% colMeans()


expressed <- tpm %>% 
  select(-target_id) %>% 
  summarise_all(funs(count = sum(.>=0), 
                     countGT0 = sum(. >0), 
                     countGTeq0.04 = sum(. >=0.04), 
                     countGTeq0.06 = sum(. >=0.06),
                     countGTeq0.1 = sum(. >=0.1),
                     countGTeq0.2 = sum(. >=0.2),
                     countGTeq0.3 = sum(. >=0.3),
                     countGTeq0.4 = sum(. >=0.4),
                     countGTeq0.5 = sum(. >=0.5),
                     countGTeq0.6 = sum(. >=0.6),
                     countGTeq0.7 = sum(. >=0.7),
                     countGTeq0.8 = sum(. >=0.8),
                     countGTeq0.9 = sum(. >=0.9),
                     countGTeq1 = sum(. >=1)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=196),
	rep(LETTERS[1:2],times=14),
	rep(LETTERS[1:3],times=42),
	rep(LETTERS[1:2],times=14),
	rep(LETTERS[1:1],times=14) )) 




expressedPlot <- ggplot() + 
  geom_line(data=subset(expressed,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_line(data=subset(expressed,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep)  ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressed,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) +
  geom_hline(yintercept=37223,color="grey50")

ggsave(expressedPlot, file="out-01.expressedPlot.jpg", width=21, height=6, units="in")

expressedPlotPoints <- ggplot() + 
  geom_point(data=subset(expressed,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_point(data=subset(expressed,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, )  ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressed,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic,) ) +
  geom_hline(yintercept=37223,color="grey50")

ggsave(expressedPlotPoints, file="out-01.expressedPlotPoints.jpg", width=21, height=6, units="in")


expressedTable <- expressed %>%
  mutate(sample=paste0(DPA,"_",rep)) %>%
  select(sample,statistic,V1) %>%
  spread(sample,V1)

write.table(file="out-01.expressedTable.tsv", expressedTable, quote=F, row.names=F, sep="\t")

sampleNames <- expressed %>% 
    mutate(sample = paste0(DPA,"_",gsub("_..DPA.*","",row.names(.)),".",rep)) %>%
    mutate(sample = gsub("AD1_TM1_","",sample)) %>%
    select(-statistic) %>%
    select(-DPA) %>% 
    select(-V1) %>% 
    select(-rep) %>%
    `rownames<-`( NULL ) %>%
    distinct()

#####################################################
##### generate similar information for only TFs #####
#####################################################

TFdata <- fread("TF.list", header=T)
TF <- TFdata$homoeolog

TFexpression <- subset(tpm, target_id %in% TF)


TFwithZeroes <- TFexpression %>% 
  select(-target_id) %>% 
  summarise_all(funs(min = min(.), 
                               q25 = quantile(., 0.25), 
                               median = median(.), 
                               medianGT0 = median(.[.>0]), 
                               medianq5.95 = median(.[.> quantile(., 0.05) & .< quantile(., 0.95)]), 
                               q75 = quantile(., 0.75), 
                               max = max(.), 
                               mean = mean(.), 
                               meanGT0 = mean(.[.>0]),
                               meanq5.95 = mean(.[.> quantile(., 0.05) & .< quantile(., 0.95)]), 
                               sd = sd(.), 
                               q5 = quantile(., 0.05), 
                               q95 = quantile(., 0.95)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=182),
	rep(LETTERS[1:2],times=13),
	rep(LETTERS[1:3],times=39),
	rep(LETTERS[1:2],times=13),
	rep(LETTERS[1:1],times=13) ))

TFsummaryPlot <- ggplot() + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "mean"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "median"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dashed" ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "q5"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "q25"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "q75"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dotted" ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "q95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "meanGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "medianGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "meanq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroes,statistic %in% "medianq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  scale_color_brewer(palette = "Dark2")

ggsave(TFsummaryPlot, file="out-01.TFsummaryPlot.jpg", width=10, height=5, units="in")

TFsummaryTable <- TFwithZeroes %>% 
  mutate(sample=paste0(DPA,"_",rep)) %>%
  select(sample,statistic,V1) %>%
  spread(sample,V1)

write.table(file="out-01.TFsummaryTable.tsv", TFsummaryTable, quote=F, row.names=F, sep="\t")



TFexpression %>% 
  select(-target_id) %>% 
  summarise_all(funs(q5GT0 = quantile(.[.>0],0.05)), 
	            	 na.rm = TRUE) %>% 
     t() %>% colMeans()


TFexpressed <- TFexpression %>% 
  select(-target_id) %>% 
  summarise_all(funs(count = sum(.>=0), 
                     countGT0 = sum(. >0), 
                     countGTeq0.025 = sum(. >=0.025),
                     countGTeq0.05 = sum(. >=0.05),
                     countGTeq0.075 = sum(. >=0.075),
                     countGTeq0.1 = sum(. >=0.1),
                     countGTeq0.2 = sum(. >=0.2),
                     countGTeq0.3 = sum(. >=0.3),
                     countGTeq0.4 = sum(. >=0.4),
                     countGTeq0.5 = sum(. >=0.5),
                     countGTeq0.6 = sum(. >=0.6),
                     countGTeq0.7 = sum(. >=0.7),
                     countGTeq0.8 = sum(. >=0.8),
                     countGTeq0.9 = sum(. >=0.9),
                     countGTeq1 = sum(. >=1)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=210),
	rep(LETTERS[1:2],times=15),
	rep(LETTERS[1:3],times=45),
	rep(LETTERS[1:2],times=15),
	rep(LETTERS[1:1],times=15) )) 




TFexpressedPlot <- ggplot() + 
  geom_line(data=subset(TFexpressed,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep)  ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressed,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) +
  geom_hline(yintercept=4932,color="grey50") + 
  geom_hline(yintercept=2466,color="grey50")

ggsave(TFexpressedPlot, file="out-01.TFexpressedPlot.jpg", width=21, height=6, units="in")

TFexpressedPlotPoints <- ggplot() + 
  geom_point(data=subset(TFexpressed,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, )  ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressed,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic,) ) +
  geom_hline(yintercept=4932,color="grey50") + 
  geom_hline(yintercept=2466,color="grey50")

ggsave(TFexpressedPlotPoints, file="out-01.TFexpressedPlotPoints.jpg", width=21, height=6, units="in")


TFexpressedTable <- TFexpressed %>%
  mutate(sample=paste0(DPA,"_",rep)) %>%
  select(sample,statistic,V1) %>%
  spread(sample,V1)

write.table(file="out-01.TFexpressedTable.tsv", TFexpressedTable, quote=F, row.names=F, sep="\t")

TFsampleNames <- TFexpressed %>% 
    mutate(sample = paste0(DPA,"_",gsub("_..DPA.*","",row.names(.)),".",rep)) %>%
    mutate(sample = gsub("AD1_TM1_","",sample)) %>%
    select(-statistic) %>%
    select(-DPA) %>% 
    select(-V1) %>% 
    select(-rep) %>%
    `rownames<-`( NULL ) %>%
    distinct()


################################################################
################################################################
#####       Sample 14 DPA, plant 11 appears weird.         #####
#####            I think we should remove it and 25 DPA    #####
################################################################
################################################################


withZeroesFiltered <- withZeroes %>% 
  filter(!(DPA=="14DPA" & rep=="A")) %>%
  filter(!(DPA=="25DPA"))

summaryPlotFil <- ggplot() + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "mean"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "median"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dashed" ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "q5"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "q25"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "q75"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dotted" ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "q95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "meanGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "medianGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "meanq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(withZeroesFiltered,statistic %in% "medianq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  scale_color_brewer(palette = "Dark2")

ggsave(summaryPlotFil, file="out-01.summaryPlotFiltered25DPA.jpg", width=6, height=3, units="in")


expressedFil <- tpm %>% 
  select(-c(target_id,AD1_TM1_plant11_14DPA,AD1_TM1_plant6_25DPA)) %>% 
  summarise_all(funs(count = sum(.>=0), 
                     countGT0 = sum(. >0), 
                     countGTeq0.04 = sum(. >=0.04), 
                     countGTeq0.06 = sum(. >=0.06),
                     countGTeq0.1 = sum(. >=0.1),
                     countGTeq0.2 = sum(. >=0.2),
                     countGTeq0.3 = sum(. >=0.3),
                     countGTeq0.4 = sum(. >=0.4),
                     countGTeq0.5 = sum(. >=0.5),
                     countGTeq0.6 = sum(. >=0.6),
                     countGTeq0.7 = sum(. >=0.7),
                     countGTeq0.8 = sum(. >=0.8),
                     countGTeq0.9 = sum(. >=0.9),
                     countGTeq1 = sum(. >=1)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=112),
	rep(LETTERS[2:3],times=14),
	rep(LETTERS[1:3],times=70),
	rep(LETTERS[1:2],times=14),
	rep(LETTERS[1:3],times=42),
	rep(LETTERS[1:2],times=14) )) 




expressedPlotFil <- ggplot() + 
  geom_line(data=subset(expressedFil,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep)  ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(expressedFil,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) +
  geom_hline(yintercept=37223,color="grey50")

ggsave(expressedPlotFil, file="out-01.expressedPlotFiltered25DPA.jpg", width=21, height=6, units="in")

expressedPlotPointsFil <- ggplot() + 
  geom_point(data=subset(expressedFil,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, )  ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(expressedFil,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic,) ) +
  geom_hline(yintercept=37223,color="grey50")

ggsave(expressedPlotPointsFil, file="out-01.expressedPlotPointsFiltered25DPA.jpg", width=21, height=6, units="in")

###############################################################################
##### generate similar information for only TFs without 14/25 DPA samples #####
###############################################################################

TFwithZeroesFiltered <- TFwithZeroes %>% 
  filter(!(DPA=="14DPA" & rep=="A")) %>%
  filter(!(DPA=="25DPA")) 


TFsummaryPlotFil <- ggplot() + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "mean"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "median"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dashed" ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "q5"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "q25"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "q75"), aes(x=DPA,y=V1,group=rep,color=statistic), linetype="dotted" ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "q95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "meanGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "medianGTO"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "meanq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  geom_line(data=subset(TFwithZeroesFiltered,statistic %in% "medianq5.95"), aes(x=DPA,y=V1,group=rep,color=statistic) ) + 
  scale_color_brewer(palette = "Dark2")

ggsave(TFsummaryPlotFil, file="out-01.TFsummaryPlotFiltered25DPA.jpg", width=10, height=5, units="in")

TFexpressedFil <- TFexpression %>% 
  select(-c(target_id,AD1_TM1_plant11_14DPA,AD1_TM1_plant6_25DPA)) %>% 
  summarise_all(funs(count = sum(.>=0), 
                     countGT0 = sum(. >0), 
                     countGTeq0.025 = sum(. >=0.025),
                     countGTeq0.05 = sum(. >=0.05),
                     countGTeq0.075 = sum(. >=0.075),
                     countGTeq0.1 = sum(. >=0.1),
                     countGTeq0.2 = sum(. >=0.2),
                     countGTeq0.3 = sum(. >=0.3),
                     countGTeq0.4 = sum(. >=0.4),
                     countGTeq0.5 = sum(. >=0.5),
                     countGTeq0.6 = sum(. >=0.6),
                     countGTeq0.7 = sum(. >=0.7),
                     countGTeq0.8 = sum(. >=0.8),
                     countGTeq0.9 = sum(. >=0.9),
                     countGTeq1 = sum(. >=1)), 
					 na.rm = TRUE) %>%  
  t() %>% 
  as.data.frame() %>%
  mutate(statistic = gsub("AD1.*DPA_","",row.names(.))) %>%
  mutate(statistic = gsub("._","",statistic)) %>%
  mutate(DPA = gsub("AD1.*_plant.._","",row.names(.))) %>%
  mutate(DPA = gsub("AD1.*_plant._","",DPA)) %>%
  mutate(DPA = gsub("_.*","",DPA)) %>%
  arrange(DPA) %>%
  mutate(rep=c(rep(LETTERS[1:3],times=120),
	rep(LETTERS[2:3],times=15),
	rep(LETTERS[1:3],times=75),
	rep(LETTERS[1:2],times=15),
	rep(LETTERS[1:3],times=45),
	rep(LETTERS[1:2],times=15) )) 



TFexpressedPlotFil <- ggplot() + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep)  ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) + 
  geom_line(data=subset(TFexpressedFil,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic, linetype=rep) ) +
  geom_hline(yintercept=4932,color="grey50") + 
  geom_hline(yintercept=2466,color="grey50")

ggsave(TFexpressedPlot, file="out-01.TFexpressedPlot25DPA.jpg", width=21, height=6, units="in")

TFexpressedPlotPointsFil <- ggplot() + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "count"), aes(x=DPA,y=V1,group=rep,color=statistic),color="black" ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGT0"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.04"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.06"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.1"), aes(x=DPA,y=V1,group=rep,color=statistic, )  ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.2"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.3"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.4"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.5"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.6"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.7"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.8"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq0.9"), aes(x=DPA,y=V1,group=rep,color=statistic, ) ) + 
  geom_point(data=subset(TFexpressedFil,statistic %in% "countGTeq1"), aes(x=DPA,y=V1,group=rep,color=statistic,) ) +
  geom_hline(yintercept=4932,color="grey50") + 
  geom_hline(yintercept=2466,color="grey50")

ggsave(TFexpressedPlotPoints, file="out-01.TFexpressedPlotPoints25DPA.jpg", width=21, height=6, units="in")

