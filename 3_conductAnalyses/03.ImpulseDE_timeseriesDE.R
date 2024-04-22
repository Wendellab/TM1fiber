library(tidyverse)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(ggplotify)
library(data.table)
library(ggrepel)
library(ggforce)


# turn off scientific notation and set seed
options(scipen = 999)
set.seed(8675309)
setDTthreads(12) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
max.overlaps=Inf


c25 <- c(
   "dodgerblue2", "#E31A1C", # red
   "green4",
   "#6A3D9A", # purple
   "#FF7F00", # orange
   "black", "orange4",
   "navy", "deeppink4", 
   "seagreen",
   "purple4", 
   "sienna2", 
   "gray30", 
   "maroon", "deeppink1", "blue1", "steelblue4",
   "darkturquoise", "green1", "yellow4", "yellow3",
   "darkorange4", "brown"
)


######### ######### ######### ######### ######### ######### ######### 
######### read in files and prepare datasets #### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 

setwd("W:/corrinne/fiberDan/RNAseq/pseudoAD1/DEanalysis")



##### ##### ##### ##### ##### ##### 
##### read in files ### ##### #####
##### ##### ##### ##### ##### ##### 

countdata <- fread("count.Purdue.pseudoAD1.filtered.25DPA.tsv", sep="\t")
metadata <- fread("metadata.Purdue.pseudoAD1.filtered.25DPA.tsv", sep="\t")


##### ##### ##### ##### ##### 
##### DGE by ImpulseDE2 ##### 
##### ##### ##### ##### ##### 

library(devtools)
install_github("YosefLab/ImpulseDE2")

library(ImpulseDE2)

# I renamed columns to match what they had in the recommendations
# I also had to rename the timepoints numerically

dfAnnot <- metadata %>%
	dplyr::rename(Sample = id, Condition = species, Time = DPA, Batch = rep) %>%
	dplyr::select(Sample,Condition,Time,Batch) %>%
	mutate(Condition = gsub("AD1","case",Condition)) %>%
	mutate(Batch = "B_NULL") %>%
	mutate(Time = gsub("06DPA",as.numeric("1"),Time)) %>% 
	mutate(Time = gsub("07DPA",as.numeric("2"),Time)) %>% 
	mutate(Time = gsub("08DPA",as.numeric("3"),Time)) %>% 
	mutate(Time = gsub("09DPA",as.numeric("4"),Time)) %>% 
	mutate(Time = gsub("10DPA",as.numeric("5"),Time)) %>% 
	mutate(Time = gsub("11DPA",as.numeric("6"),Time)) %>% 
	mutate(Time = gsub("12DPA",as.numeric("7"),Time)) %>% 
	mutate(Time = gsub("13DPA",as.numeric("8"),Time)) %>% 
	mutate(Time = gsub("14DPA",as.numeric("9"),Time)) %>% 
	mutate(Time = gsub("15DPA",as.numeric("10"),Time)) %>% 
	mutate(Time = gsub("16DPA",as.numeric("11"),Time)) %>% 
	mutate(Time = gsub("17DPA",as.numeric("12"),Time)) %>% 
	mutate(Time = gsub("18DPA",as.numeric("13"),Time)) %>% 
	mutate(Time = gsub("19DPA",as.numeric("14"),Time)) %>% 
	mutate(Time = gsub("20DPA",as.numeric("15"),Time)) %>% 
	mutate(Time = gsub("21DPA",as.numeric("16"),Time)) %>% 
	mutate(Time = gsub("22DPA",as.numeric("17"),Time)) %>% 
	mutate(Time = gsub("23DPA",as.numeric("18"),Time)) %>% 
	mutate(Time = gsub("24DPA",as.numeric("19"),Time)) %>% 
	mutate_at('Time',as.numeric) %>% 
	arrange(Time)

# I also had to take the count df and adjust it to their parameters
# genes become rownames
# 

matCount <- countdata %>%
	column_to_rownames(var="target_id") %>%
	dplyr::select(dfAnnot$Sample) %>% 
	as.matrix
	
IDE2 <- runImpulseDE2(
	matCountData = matCount,
	dfAnnotation = dfAnnot,
	boolCaseCtrl = FALSE,
	boolIdentifyTransients = TRUE,
	vecConfounders = NULL,
	scaNProc = 4)
	
save(IDE2, file = "03.R-ImpulseDE2.RData")


IDE2results <- as.data.frame(IDE2$dfImpulseDE2Results)
write.table(IDE2results,file="out-03.IDE2results.filtered.25DPA.tsv",quote=F, sep="\t", row.names=F)




##### ##### ##### ##### ##### ##### #####
##### read in CesA genes             #### 
##### these are my genes of interest ####
##### ##### ##### ##### ##### ##### #####

PStable <- read.table("PCWSCW.tsv", header =T)
PCWtable <- PStable %>%
	filter(stage =="PCW")

SCWtable <- PStable %>%
	filter(stage =="SCW")

CesA <- as.vector(PStable$homoeolog[PStable$homoeolog %in% rownames(IDE2@matCountDataProc)])
PCW <- as.vector(PCWtable$homoeolog[PCWtable$homoeolog %in% rownames(IDE2@matCountDataProc)])
SCW <- as.vector(SCWtable$homoeolog[SCWtable$homoeolog %in% rownames(IDE2@matCountDataProc)])


### impulse profiles for CesA genes
CesAImpulse <- plotGenes(objectImpulseDE2 = IDE2, vecGeneIDs=CesA, boolCaseCtrl=F, dirOut='./', strFileName="out-03.Impulse_CesA.pdf",boolMultiplePlotsPerPage=T)
PCWImpulse <- plotGenes(objectImpulseDE2 = IDE2, vecGeneIDs=PCW, boolCaseCtrl=F, dirOut='./', strFileName="out-03.Impulse_PCW.pdf",boolMultiplePlotsPerPage=T)
SCWImpulse <- plotGenes(objectImpulseDE2 = IDE2, vecGeneIDs=SCW, boolCaseCtrl=F, dirOut='./', strFileName="out-03.Impulse_SCW.pdf",boolMultiplePlotsPerPage=T)




##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### potential visualization and gene lists ## ##### 
##### ##### ##### ##### ##### ##### ##### ##### ##### 

library(ComplexHeatmap)

# makes the nice gradient picture
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = IDE2,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = 0.01)
draw(lsHeatmaps$complexHeatmapRaw) 

jpeg(filename = "out-03.complexHeatmapFit.jpg",
	width = 8,
	height = 8,
	res = 300,
	units = "in",
	pointsize = 12)

draw(lsHeatmaps$complexHeatmapFit) 

dev.off()


# gives you what is transient/transition up/down
lsHeatmaps$lsvecGeneGroups

# gives you the number in each category
length(lsHeatmaps$lsvecGeneGroups$transient_up)
#[1] 3402
length(lsHeatmaps$lsvecGeneGroups$transient_down)
#[1] 1871
length(lsHeatmaps$lsvecGeneGroups$transition_up)
#[1] 19076
length(lsHeatmaps$lsvecGeneGroups$transition_down)
#[1] 14491

# writes gene lists for each category, but you can search for genes of interest
# before you do any writing
fwrite(list(lsHeatmaps$lsvecGeneGroups$transient_up), file="out-03.transient_up.list", sep="\n",col.names=F)
fwrite(list(lsHeatmaps$lsvecGeneGroups$transient_down), file="out-03.transient_down.list", sep="\n",col.names=F)
fwrite(list(lsHeatmaps$lsvecGeneGroups$transition_up), file="out-03.transition_up.list", sep="\n",col.names=F)
fwrite(list(lsHeatmaps$lsvecGeneGroups$transition_down), file="out-03.transition_down.list", sep="\n",col.names=F)


Iup <- scan("out-03.transient_up.list", what="character")
Idown <- scan("out-03.transient_down.list", what="character")
Tup <- scan("out-03.transition_up.list", what="character")
Tdown <- scan("out-03.transition_down.list", what="character")


### I saved everything, but then I read it back in to search it for 
### genes of interest, here transcription factors

TFs <- fread("GoraiTFs.txt",  header=T)
TFs <- bind_rows(TFs %>% mutate(Gene_ID = gsub("$",".1.A",Gene_ID)),TFs %>% mutate(Gene_ID = gsub("$",".1.D",Gene_ID)))

IupTF <- TFs %>% filter(Gene_ID %in% Iup)
IdownTF <- TFs %>% filter(Gene_ID %in% Idown)
TupTF <- TFs %>% filter(Gene_ID %in% Tup)
TdownTF <- TFs %>% filter(Gene_ID %in% Tdown)


fwrite(IupTF, file="out-03.transient_up.TF.list", sep="\n",col.names=F)
fwrite(IdownTF, file="out-03.transient_down.TF.list", sep="\n",col.names=F)
fwrite(TupTF, file="out-03.transition_up.TF.list", sep="\n",col.names=F)
fwrite(TdownTF, file="out-03.transition_down.TF.list", sep="\n",col.names=F)


nrow(IupTF)
nrow(IupTF)/length(Iup)
# 343
# 0.100823


nrow(IdownTF)
nrow(IdownTF)/length(Idown)
# 334
# 0.1785142


nrow(TupTF)
nrow(TupTF)/length(Tup)
# 3868
# 0.2027679


nrow(TdownTF)
nrow(TdownTF)/length(Tdown)
# 1540
# 0.1062729

# proportion of TF that are impulse different from constant for each direction
prop.test(343,3402,p=0.203)
prop.test(334,1871,p=0.106)




#########################
##### GO enrichment #####
#########################

load('D5annotation.Rdata')
# "annotation221" "annot"

load('cottonGOenrich.RData')
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"  "level1_BP_terms"
# "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" 
# "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" 
# "level5_MF_terms"


### Annotate groups with topGO ###

universe <- row.names(matCount)

# for each group containing genes of interest
GOresults <- data.frame(GO.ID=character(),
  Term=character(), Annotated=integer(),
  Significant=integer(), Expected=numeric(),
  result1=character(), qval.bh=numeric(),
  ontology=character(), category=character(), 
  stringsAsFactors=F)

# modify geneID2GO
a <- geneID2GO
d <- geneID2GO
names(a) <- paste0(names(a),".1.A")
names(d) <- paste0(names(d),".1.D")
geneID2GOp <- c(a,d)

dir.create("topGO")

source("function_topGo.R")

do_topGO(Iup,"impulseUP")
do_topGO(Idown,"impulseDOWN")
do_topGO(Tup,"transitionUP")
do_topGO(Tdown,"transitionDOWN")

fwrite(na.omit(GOresults), file="out-03.impulse_topGO_enrichment_filtered.txt", sep="\t", row.names=FALSE)


################################
### visualize revigo results ###
################################

source("Revigo_ImpulseBuildSource.R")


IupBP<- ggplot(IupBP.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Impulse up, Biological Process, ",length(Iup)," genes")) +
  geom_treemap(layout="squarified") + 
  geom_treemap_text(place = "centre",size = 12,reflow=T) 

IdownBP<- ggplot(IdownBP.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Impulse down, Biological Process, ",length(Idown)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)

TupBP<- ggplot(TupBP.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Transition up, Biological Process, ",length(Tup)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)

TdownBP<- ggplot(TdownBP.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Transition down, Biological Process, ",length(Tdown)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)


IupMF<- ggplot(IupMF.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Impulse up, Molecular Function, ",length(Iup)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)

IdownMF<- ggplot(IdownMF.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Impulse down, Molecular Function, ",length(Idown)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)

TupMF<- ggplot(TupMF.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Transition up, Molecular Function, ",length(Tup)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)

TdownMF<- ggplot(TdownMF.stuff, aes(area=value, fill=representative,label=description,subgroup=representative)) +
  ggtitle(paste0("     Transition down, Molecular Function, ",length(Tdown)," genes")) +
  geom_treemap(layout="squarified")+ 
  geom_treemap_text(place = "centre",size = 12,reflow=T)





IupMFS <- ggplot( data = IupMF.one.data ) + 
  ggtitle(paste0("     Impulse up, Molecular Function, ",length(Iup)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(IupMF.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = IupMF.one.data [ IupMF.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(IupMF.one.data$plot_X)-(max(IupMF.one.data$plot_X) - min(IupMF.one.data$plot_X))/10,max(IupMF.one.data$plot_X)+(max(IupMF.one.data$plot_X) - min(IupMF.one.data$plot_X))/10) + 
  ylim(min(IupMF.one.data$plot_Y)-(max(IupMF.one.data$plot_Y) - min(IupMF.one.data$plot_Y))/10,max(IupMF.one.data$plot_Y)+(max(IupMF.one.data$plot_Y) - min(IupMF.one.data$plot_Y))/10)



IdownMFS <- ggplot( data = IdownMF.one.data ) + 
  ggtitle(paste0("     Impulse down, Molecular Function", length(Idown)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(IdownMF.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = IdownMF.one.data [ IdownMF.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(IdownMF.one.data$plot_X)-(max(IdownMF.one.data$plot_X) - min(IdownMF.one.data$plot_X))/10,max(IdownMF.one.data$plot_X)+(max(IdownMF.one.data$plot_X) - min(IdownMF.one.data$plot_X))/10) + 
  ylim(min(IdownMF.one.data$plot_Y)-(max(IdownMF.one.data$plot_Y) - min(IdownMF.one.data$plot_Y))/10,max(IdownMF.one.data$plot_Y)+(max(IdownMF.one.data$plot_Y) - min(IdownMF.one.data$plot_Y))/10)


TupMFS <- ggplot( data = TupMF.one.data ) + 
  ggtitle(paste0("     Transition up, Molecular Function, ",length(Tup)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(TupMF.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = TupMF.one.data [ TupMF.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(TupMF.one.data$plot_X)-(max(TupMF.one.data$plot_X) - min(TupMF.one.data$plot_X))/10,max(TupMF.one.data$plot_X)+(max(TupMF.one.data$plot_X) - min(TupMF.one.data$plot_X))/10) + 
  ylim(min(TupMF.one.data$plot_Y)-(max(TupMF.one.data$plot_Y) - min(TupMF.one.data$plot_Y))/10,max(TupMF.one.data$plot_Y)+(max(TupMF.one.data$plot_Y) - min(TupMF.one.data$plot_Y))/10)



TdownMFS <- ggplot( data = TdownMF.one.data ) + 
  ggtitle(paste0("     Transition down, Molecular Function, ",length(Tdown)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(TdownMF.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = TdownMF.one.data [ TdownMF.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(TdownMF.one.data$plot_X)-(max(TdownMF.one.data$plot_X) - min(TdownMF.one.data$plot_X))/10,max(TdownMF.one.data$plot_X)+(max(TdownMF.one.data$plot_X) - min(TdownMF.one.data$plot_X))/10) + 
  ylim(min(TdownMF.one.data$plot_Y)-(max(TdownMF.one.data$plot_Y) - min(TdownMF.one.data$plot_Y))/10,max(TdownMF.one.data$plot_Y)+(max(TdownMF.one.data$plot_Y) - min(TdownMF.one.data$plot_Y))/10)





IupBPS <- ggplot( data = IupBP.one.data ) + 
  ggtitle(paste0("     Impulse up, Biological Process, ",length(Iup)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(IupBP.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = IupBP.one.data [ IupBP.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(IupBP.one.data$plot_X)-(max(IupBP.one.data$plot_X) - min(IupBP.one.data$plot_X))/10,max(IupBP.one.data$plot_X)+(max(IupBP.one.data$plot_X) - min(IupBP.one.data$plot_X))/10) + 
  ylim(min(IupBP.one.data$plot_Y)-(max(IupBP.one.data$plot_Y) - min(IupBP.one.data$plot_Y))/10,max(IupBP.one.data$plot_Y)+(max(IupBP.one.data$plot_Y) - min(IupBP.one.data$plot_Y))/10)



IdownBPS <- ggplot( data = IdownBP.one.data ) + 
  ggtitle(paste0("     Impulse down, Biological Process, ",length(Idown)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(IdownBP.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = IdownBP.one.data [ IdownBP.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(IdownBP.one.data$plot_X)-(max(IdownBP.one.data$plot_X) - min(IdownBP.one.data$plot_X))/10,max(IdownBP.one.data$plot_X)+(max(IdownBP.one.data$plot_X) - min(IdownBP.one.data$plot_X))/10) + 
  ylim(min(IdownBP.one.data$plot_Y)-(max(IdownBP.one.data$plot_Y) - min(IdownBP.one.data$plot_Y))/10,max(IdownBP.one.data$plot_Y)+(max(IdownBP.one.data$plot_Y) - min(IdownBP.one.data$plot_Y))/10)


TupBPS <- ggplot( data = TupBP.one.data ) + 
  ggtitle(paste0("     Transition up, Biological Process, ",length(Tup)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(TupBP.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = TupBP.one.data [ TupBP.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(TupBP.one.data$plot_X)-(max(TupBP.one.data$plot_X) - min(TupBP.one.data$plot_X))/10,max(TupBP.one.data$plot_X)+(max(TupBP.one.data$plot_X) - min(TupBP.one.data$plot_X))/10) + 
  ylim(min(TupBP.one.data$plot_Y)-(max(TupBP.one.data$plot_Y) - min(TupBP.one.data$plot_Y))/10,max(TupBP.one.data$plot_Y)+(max(TupBP.one.data$plot_Y) - min(TupBP.one.data$plot_Y))/10)



TdownBPS <- ggplot( data = TdownBP.one.data ) + 
  ggtitle(paste0("     Transition down, Biological Process, ",length(Tdown)," genes")) +
  geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) +
  scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(TdownBP.one.data$value), 0) ) + 
  geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + 
  scale_size( range=c(5, 30)) + 
  theme_bw() + 
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  geom_text( data = TdownBP.one.data [ TdownBP.one.data$dispensability < 0.15, ], aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 ) +
  labs (y = "semantic space x", x = "semantic space y") +
  theme(legend.key = element_blank()) +
  xlim(min(TdownBP.one.data$plot_X)-(max(TdownBP.one.data$plot_X) - min(TdownBP.one.data$plot_X))/10,max(TdownBP.one.data$plot_X)+(max(TdownBP.one.data$plot_X) - min(TdownBP.one.data$plot_X))/10) + 
  ylim(min(TdownBP.one.data$plot_Y)-(max(TdownBP.one.data$plot_Y) - min(TdownBP.one.data$plot_Y))/10,max(TdownBP.one.data$plot_Y)+(max(TdownBP.one.data$plot_Y) - min(TdownBP.one.data$plot_Y))/10)


ggsave(filename="topGO/IupBP_wordmap.jpg", device="jpeg", plot=IupBP, width=12, height=12)
ggsave(filename="topGO/IupMF_wordmap.jpg", device="jpeg", plot=IupMF, width=12, height=12)
ggsave(filename="topGO/IdownBP_wordmap.jpg", device="jpeg", plot=IdownBP, width=12, height=12)
ggsave(filename="topGO/IdownMF_wordmap.jpg", device="jpeg", plot=IdownMF, width=12, height=12)
ggsave(filename="topGO/TupBP_wordmap.jpg", device="jpeg", plot=TupBP, width=12, height=12)
ggsave(filename="topGO/TupMF_wordmap.jpg", device="jpeg", plot=TupMF, width=12, height=12)
ggsave(filename="topGO/TdownBP_wordmap.jpg", device="jpeg", plot=TdownBP, width=12, height=12)
ggsave(filename="topGO/TdownMF_wordmap.jpg", device="jpeg", plot=TdownMF, width=12, height=12)

ggsave(filename="topGO/IupBP_Scatter.jpg", device="jpeg", plot=IupBPS, width=12, height=12)
ggsave(filename="topGO/IupMF_Scatter.jpg", device="jpeg", plot=IupMFS, width=12, height=12)
ggsave(filename="topGO/IdownBP_Scatter.jpg", device="jpeg", plot=IdownBPS, width=12, height=12)
ggsave(filename="topGO/IdownMF_Scatter.jpg", device="jpeg", plot=IdownMFS, width=12, height=12)
ggsave(filename="topGO/TupBP_Scatter.jpg", device="jpeg", plot=TupBPS, width=12, height=12)
ggsave(filename="topGO/TupMF_Scatter.jpg", device="jpeg", plot=TupMFS, width=12, height=12)
ggsave(filename="topGO/TdownBP_Scatter.jpg", device="jpeg", plot=TdownBPS, width=12, height=12)
ggsave(filename="topGO/TdownMF_Scatter.jpg", device="jpeg", plot=TdownMFS, width=12, height=12)

BPplot <- plot_grid(IupBP,IdownBP,TupBP,TdownBP,labels='AUTO', label_size=20)
MFplot <- plot_grid(IupMF,IdownMF,TupMF,TdownMF,labels='AUTO', label_size=20)
BPSplot <- plot_grid(IupBPS,IdownBPS,TupBPS,TdownBPS,labels='AUTO', label_size=20)
MFSplot <- plot_grid(IupMFS,IdownMFS,TupMFS,TdownMFS,labels='AUTO', label_size=20)


ggsave(filename="BPwordmap.jpg", device="jpeg", plot=BPplot, width=25, height=20)
ggsave(filename="MFwordmap.jpg", device="jpeg", plot=MFplot, width=25, height=20)
ggsave(filename="BPScatter.jpg", device="jpeg", plot=BPSplot, width=22, height=20)
ggsave(filename="MFScatter.jpg", device="jpeg", plot=MFSplot, width=22, height=20)


