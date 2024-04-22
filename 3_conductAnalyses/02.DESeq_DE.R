library(tidyverse)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(ggplotify)
library(data.table)
library(ggrepel)
library(DEGreport)
library(ggforce)

##### ##### ##### ##### ##### ##### 
##### establish functions ### ##### 
##### ##### ##### ##### ##### ##### 

# turn off scientific notation and set seed
options(scipen = 999)
set.seed(8675309)
setDTthreads(12) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
max.overlaps=Inf

#### function to merge tsv files from kallisto, report counts #### 
source("function_KallistoMerge.R")

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

countdata <- fread("count.Purdue.pseudoAD1.tsv", sep="\t")

tpmdata <- fread("tpm.Purdue.pseudoAD1.tsv", sep="\t")
metadata <- fread("metadata.Purdue.pseudoAD1.tsv", sep="\t")

# bad sample from 01.ExpressionSummary.R is AD1_TM1_plant11_14DPA_tpm
tpmMean <- tpmdata %>% 
  select(-AD1_TM1_plant11_14DPA_tpm) %>% 
  gather(sample, value, -target_id) %>% 
  mutate(sample = gsub("plant.._","",sample)) %>% 
  mutate(sample = gsub("plant._","",sample)) %>% 
  mutate(sample = gsub("_._","_",sample)) %>%
  group_by(target_id,sample) %>% 
  summarize_at(vars(value), list(value=mean)) %>%
  spread(sample,value)

write.table(tpmMean, file="tpmMean.Purdue.pseudoAD1.clean.tsv", quote=F, sep="\t", row.names=F)





##### ##### ##### ##### ##### ##### 
##### Visualize all samples # ##### 
##### ##### ##### ##### ##### ##### 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~DPA, tidy = TRUE)
dds <- dds[rowSums(counts(dds))/nrow(metadata) >= 1,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)


### make PCA/pheatmap images ###
vsd <- vst(dds, blind=FALSE)

customPCA <- plotPCA(vsd, intgroup=c("DPA"), returnData=T, ntop=Inf)
percentVar <- round(100 * attr(customPCA, "percentVar"))
reviewPCA <- ggplot(customPCA, aes(PC1, PC2, color=DPA)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()

ggsave("out-02.reviewPCA.pseudoAD1.jpg", plot=reviewPCA, width=18, height=12, units="in")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$DPA, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$DPA, sep="-")
c <- as.grob(pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,silent=T))

myplot <- plot_grid(reviewPCA,c)

ggsave("out-02.PCA_pheatmap.pseudoAD1.jpg",plot=myplot, width=20, height=10, units="in")

### umap visualization ###
library(umap)

fiberFeatures <- t(as.data.frame(assay(vsd))) 
fiber.umap <- umap(fiberFeatures)
fiber.layout <- data.frame(fiber.umap[["layout"]]) 
fiber.layout$samples <- row.names(fiber.layout)
fiber.layout$DPA <- gsub("AD1_TM1_.*_","",fiber.layout$samples)

umapHplot <- ggplot(fiber.layout, aes(X1, X2, color=DPA)) +
  geom_point(size=3) +
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()

ggsave("out-02.umapPlot.homoeolog.pseudoAD1.jpg",plot=umapHplot, width=10, height=12, units="in")







##### ##### ##### ##### ##### ##### ##### ##### 
##### Visualize only filtered samples ### ##### 
##### ##### ##### ##### ##### ##### ##### #####

### remove AD1_TM1_plant11_14DPA, AD1_TM1_plant6_25DPA
remove <- c("AD1_TM1_plant11_14DPA", "AD1_TM1_plant6_25DPA")
countdata <- countdata %>% select(-c(AD1_TM1_plant11_14DPA, AD1_TM1_plant6_25DPA))
metadata <- metadata %>% filter(!id %in% remove)

write.table(countdata, file="count.Purdue.pseudoAD1.filtered.25DPA.tsv", quote=F, sep="\t", row.names=F)
write.table(metadata, file="metadata.Purdue.pseudoAD1.filtered.25DPA.tsv", quote=F, sep="\t", row.names=F)


### make new, filtered dds ###
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~DPA, tidy = TRUE)
dds <- dds[rowSums(counts(dds))/nrow(metadata) >= 1,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)


### make PCA/pheatmap images ###
vsd <- vst(dds, blind=FALSE)

customPCA <- plotPCA(vsd, intgroup=c("DPA"), returnData=T, ntop=Inf)
percentVar <- round(100 * attr(customPCA, "percentVar"))
reviewPCA <- ggplot(customPCA, aes(PC1, PC2, color=DPA)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()

ggsave("out-02.reviewPCA.pseudoAD1.filtered.25DPA.jpg", plot=reviewPCA, width=18, height=12, units="in")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$DPA, sep="-")
colnames(sampleDistMatrix) <- paste(vsd$condition, vsd$DPA, sep="-")
c <- as.grob(pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,silent=T))

myplot <- plot_grid(reviewPCA,c)

ggsave("out-02.PCA_pheatmap.pseudoAD1.filtered.25DPA.jpg",plot=myplot, width=20, height=10, units="in")

### umap visualization ###
library(umap)

fiberFeatures <- t(as.data.frame(assay(vsd))) 
fiber.umap <- umap(fiberFeatures)
fiber.layout <- data.frame(fiber.umap[["layout"]]) 
fiber.layout$samples <- row.names(fiber.layout)
fiber.layout$DPA <- gsub("AD1_TM1_.*_","",fiber.layout$samples)

umapHplot <- ggplot(fiber.layout, aes(X1, X2, color=DPA)) +
  geom_point(size=3) +
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()

ggsave("out-02.umapPlot.homoeolog.pseudoAD1.filtered.25DPA.jpg",plot=umapHplot, width=10, height=12, units="in")











##### ##### ##### ##### ##### ##### ##### #####  
#####  Differential expression #### ##### #####
##### ##### ##### ##### ##### ##### ##### ##### 


source("function_DEcontrastsDPA.R")
##### DPA contrasts #####

DPAcontrast("16","14")

DPAcontrast("07","06")
DPAcontrast("08","07")
DPAcontrast("09","08")
DPAcontrast("10","09")

# make contrasts for 11v10 through 25v24
for (i in 11:24) {
	j=i-1
	DPAcontrast(i,j)
}


write.table(contrastTable,file="out-02.Contrasts.filtered.25DPA.tsv",quote=F, sep="\t", row.names=F)


for (i in 2:ncol(SigFoldChange)) {
	up <- SigFoldChange[,c(1,i)] %>% na.omit %>% dplyr::rename(log2FC = 2) %>% filter_at(2,all_vars(. >0))
	down <- SigFoldChange[,c(1,i)] %>% na.omit %>% dplyr::rename(log2FC = 2) %>% filter_at(2,all_vars(. <0))
	write.table(up,file=paste0("DEcontrasts/out-02.",colnames(SigFoldChange[i]),".significant.up.tsv"), quote=F, sep="\t", row.names=F)
	write.table(down,file=paste0("DEcontrasts/out-02.",colnames(SigFoldChange[i]),".significant.down.tsv"), quote=F, sep="\t", row.names=F)
}

write.table(SigFoldChange,file="out-02.SigFoldChange.filtered.25DPA.tsv",quote=F, sep="\t", row.names=F)


ggpContrast <- as.data.frame(contrastTable) %>% 
	mutate(across(.cols=2:10, .fns=as.numeric)) %>% 
	pivot_longer(-contrast) %>%
	mutate(across(.cols=1:2, .fns=as.factor)) %>%
	mutate(DEgrp=name) %>% 
	mutate(name=fct_relevel(name,c("DE","DEup","DEdown","DEA","DEAup","DEAdown","DED","DEDup","DEDdown"))) %>%
  mutate(grp=gsub("DEA.*","A homoeolog",name)) %>% 
  mutate(grp=gsub("DED.*","D homoeolog",grp)) %>% 
  mutate(grp=gsub("DE.*","both homoeologs",grp)) %>%
  mutate(dir=gsub("DE.*up","up",DEgrp)) %>%
  mutate(dir=gsub("DE.*down","down",dir)) %>%
  mutate(dir=gsub("DE.*","DE",dir))

ggpContrast2 <- ggpContrast %>%
  ggplot(aes(x = contrast, y = value)) + 
  geom_line(data = transform(ggpContrast, name = NULL), aes(group = DEgrp), colour = "grey80") +
  geom_line(aes(group = name)) +
  facet_wrap(~name) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("out-02.DPAdailyContrast.filtered.25DPA.jpg", plot=ggpContrast2, width=24, height=12, units="in")

ggpContrast3 <- ggplot(ggpContrast, aes(x = contrast, y = value, color=grp)) +
  geom_line(aes(group = name, linetype = grp), size = 1) +
  facet_wrap(~dir) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 18), 
        axis.text.y = element_text(size = 18), 
        strip.text = element_text(size = 24), 
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        legend.text = element_text(size = 14),  
        legend.title = element_blank())

ggsave("out-02.DPAdailyContrast.filtered.25DPA.concise.jpg", plot=ggpContrast3, width=24, height=12, units="in")


############################################
##### summing the homoeologs into genes ####
############################################

Pcountdata <- countdata %>%
	mutate(gene=gsub(".A","",target_id)) %>%
	mutate(gene=gsub(".D","",gene)) 

Scountdata <- aggregate(. ~gene, data=Pcountdata[,-1],sum)


#### function to perform DESeq ####
Sdds <- DESeqDataSetFromMatrix(countData=Scountdata, colData=metadata, design=~DPA, tidy = TRUE)
Sdds <- Sdds[rowSums(counts(Sdds))/nrow(metadata) >= 1,]
Sdds <- estimateSizeFactors(Sdds)
Sdds <- DESeq(Sdds)


#### function to make PCA/pheatmap images ####
Svsd <- vst(Sdds, blind=FALSE)

ScustomPCA <- plotPCA(Svsd, intgroup=c("DPA"), returnData=T, ntop=Inf)
SpercentVar <- round(100 * attr(ScustomPCA, "percentVar"))
SreviewPCA <- ggplot(ScustomPCA, aes(PC1, PC2, color=DPA)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",SpercentVar[1],"% variance")) +
  ylab(paste0("PC2: ",SpercentVar[2],"% variance")) + 
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()

ggsave("out-02.reviewPCA.pseudoAD1.summed.filtered.25DPA.jpg", plot=SreviewPCA, width=18, height=12, units="in")


### summed umap ###
SfiberFeatures <- t(as.data.frame(assay(Svsd))) 
Sfiber.umap <- umap(SfiberFeatures)
Sfiber.layout <- data.frame(Sfiber.umap[["layout"]]) 
Sfiber.layout$samples <- row.names(Sfiber.layout)
Sfiber.layout$DPA <- gsub("AD1_TM1_.*_","",Sfiber.layout$samples)


umapSplot <- ggplot(Sfiber.layout, aes(X1, X2, color=DPA)) +
  geom_point(size=3) +
  coord_fixed() +
  geom_text_repel(box.padding = 0.3,aes(label=DPA)) + 
  scale_color_manual(values=c25) +
  geom_mark_ellipse()


ggsave("out-02.umapPlot.summed.pseudoAD1.filtered.25DPA.jpg",plot=umapSplot, width=10, height=12, units="in")



### make summed homoeolog contrast table ###

SDPAcontrast("07","06")
SDPAcontrast("08","07")
SDPAcontrast("09","08")
SDPAcontrast("10","09")

# make contrasts for 11v10 through 25v24
for (i in 11:24) {
	j=i-1
	SDPAcontrast(i,j)
}

write.table(ScontrastTable,file="out-02.Contrasts.summed.filtered.25DPA.tsv",quote=F, sep="\t", row.names=F)




SggpContrast <- as.data.frame(ScontrastTable) %>% 
	mutate(across(.cols=2:4, .fns=as.numeric)) %>% 
	pivot_longer(-contrast) %>%
	mutate(across(.cols=1:2, .fns=as.factor)) %>%
	mutate(DEgrp=name) %>% 
	mutate(name=fct_relevel(name,c("DE","DEup","DEdown"))) 

SggpContrast <- SggpContrast %>%
  ggplot(aes(x = contrast, y = value)) + 
  geom_line(data = transform(SggpContrast, name = NULL), aes(group = DEgrp), colour = "grey80") +
  geom_line(aes(group = name)) +
  facet_wrap(~name) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave("out-02.DPAdailyContrast.summed.filtered.25DPA.jpg", plot=SggpContrast, width=24, height=12, units="in")





