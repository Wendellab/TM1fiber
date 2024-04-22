library(stringi)

library(RColorBrewer)
library(data.table)
library(tidyverse)
library(WGCNA)

options(scipen = 999)
set.seed(8675309)
setDTthreads(18) # setting the data.table threads
options(datatable.fread.datatable=FALSE)
options(stringsAsFactors = FALSE); # do not omit
enableWGCNAThreads(nThreads=); # set WGCNA threads


#############################
## 6. network associations ##
#############################

load(file = "04.R-01-dataInput.filtered.RData")
load(file = "06.R-03-buildNetwork.RData")
load(file = "07.R-04-networkTopology.RData")


# Define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recalculate MEs with color labels
MEorder <- moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs <- orderMEs(MEorder)

sampleNames <- bind_cols(rownames(datExpr),gsub("DPA_.*","DPA",rownames(datExpr)) ) %>% 
    rename(Sample = 1, DPA = 2) %>%
    mutate(DPA = gsub("AD1.*_","",DPA)) %>% 
    arrange(DPA) %>% 
    as.data.frame()

traitData <- read.table("phenotypes/traitData.txt",sep="\t",header=T)
wallData <- traitData %>%
	dplyr::select(-c(Rep,Elemental.Growth.Rate,First.derivative.growth.rate,Microfibril..RL..orientation,Microfibril..tip..orientation)) %>%
	bind_cols(.,sampleNames$Sample) %>% 
	rename(Sample = 8) %>%
	column_to_rownames(var = "Sample") %>%
	mutate(DPA=as.numeric(sub("DPA","",DPA)))


sampleorder <- rownames(datExpr)
traitRows <- match(sampleorder, row.names(wallData))
datTraits <- wallData[traitRows,]

traitColors <- numbers2colors(datTraits, signed = F)
sampleTrees = hclust(dist(datExpr), method = "average")

pdf(file = "out-09.s1.SampleClusteringWithTraits.filtered.pdf", width = 12, height = 25);

plotDendroAndColors(sampleTrees, traitColors,
                    groupLabels = names(datTraits),
                    main = "Dendrogram and trait heatmap")


dev.off()




### correlate modules and traits ###
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Only disply correlations
sigMatrix <- round(moduleTraitCor,2)
sigMatrix[moduleTraitPvalue > 0.05] <- ""

# Name the modules by number, not color
old <- NULL
new <- NULL

for (me in 0:(Nmodules-1)) {
  old <- c(old,paste0("ME",labels2colors(me)))
  new <- c(new,paste0("ME",me))
}


MEsName <- MEs

names(MEsName) <- stri_replace_all_regex(names(MEsName),
                                  pattern= '\\b'%s+%c(old)%s+%'\\b',
                                  replacement=new,
                                  vectorize=F)



pdf(file = "out-09.s6.TraitCorrelations.filtered.pdf", width = 10, height = 10)
# Display the correlation values within a heatmap plot
par(mar=c(8,5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(moduleTraitCor) %>%
                 gsub("il..","il (",.) %>%
                 gsub("..ori", ") ori", .) %>%
                 gsub("..Mpa.", " (Mpa)", .) %>%
                 gsub("\\.", " ", .),
               yLabels = names(MEs),
               ySymbols = names(MEsName),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
		   textMatrix = sigMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



dev.off()



### define a gene significance variables ###
cellulose <- as.numeric(cor(datExpr,datTraits$Cellulose,use="p"))
celluloseColor <- numbers2colors(cellulose,signed=T)

pectin <- as.numeric(cor(datExpr,datTraits$Pectin,use="p"))
pectinColor <- numbers2colors(pectin,signed=T)

hemicellulose <- as.numeric(cor(datExpr,datTraits$Hemicellulose,use="p"))
hemicelluloseColor <- numbers2colors(hemicellulose,signed=T)

DPA <- as.numeric(cor(datExpr,datTraits$DPA,use="p"))
dpaColor <- numbers2colors(DPA,signed=T)

boll_len <- as.numeric(cor(datExpr,datTraits$Boll.Length,use="p"))
boll_lenColor <- numbers2colors(boll_len,signed=T)

fiber_len <- as.numeric(cor(datExpr,datTraits$Fiber.Length,use="p"))
fiber_lenColor <- numbers2colors(fiber_len,signed=T)

EGR <- as.numeric(cor(datExpr,datTraits$Elemental.Growth.Rate,use="p"))
EGRColor <- numbers2colors(EGR,signed=T)

FDGR <- as.numeric(cor(datExpr,datTraits$First.derivative.growth.rate,use="p"))
FDGRColor <- numbers2colors(FDGR,signed=T)

MrlO <- as.numeric(cor(datExpr,datTraits$Microfibril..RL..orientation,use="p"))
MrlOColor <- numbers2colors(MrlO,signed=T)

MtO <- as.numeric(cor(datExpr,datTraits$Microfibril..tip..orientation,use="p"))
MtOColor <- numbers2colors(MtO,signed=T)

Tgr <- as.numeric(cor(datExpr,datTraits$Turgor..Mpa,use="p"))
TgrColor <- numbers2colors(Tgr,signed=T)

pdf(file = "out-09.s6.TraitCorrelations.dendro.pdf", width = 30, height = 15)
plotDendroAndColors(net$dendrograms[[1]],
	colors=data.frame(labels2colors(net$colors),
				celluloseColor, 
				pectinColor, 
				hemicelluloseColor,
				dpaColor,
				boll_lenColor,
				fiber_lenColor,
				EGRColor,
				FDGRColor,
				MrlOColor,
				MtOColor,
				TgrColor), 
	groupLabels=c("modules","cellulose","pectin","hemicellulose","DPA", "boll length","fiber length", "exponential growth rater",
	               "first derivative growth rate", "microfibril (RL) orientation", "microfibril (tip) orientation", "turgor (Mpa)"),
	dendroLabels=FALSE, 
	hang=0.03,
	addGuide=TRUE,
	guideHang=0.05, 
	main="")

dev.off()


PStable <- read.table("PCWSCW.tsv", header =T)

exprNorm <- t(datExpr) %>% 
	as.data.frame() %>% 
	rownames_to_column("homoeolog") %>% 
	right_join(.,PStable)

CesAexpr <- exprNorm %>%
	pivot_longer(!(c(homoeolog,CesA,stage)),names_to="Sample", values_to="expression") %>% 
	left_join(.,coldata[,c("id","DPA")],by=c("Sample"="id")) %>%
	mutate(DPA = as.numeric(sub("DPA","",DPA))) %>% 
	mutate(label=ifelse(Sample=="AD1_TM1_plant10_24DPA",CesA,NA)) 

labelInfo <-
  split(CesAexpr, d$z) %>%
  lapply(function(t) {
    data.frame(
      predAtMax = loess(y ~ x, span = 0.8, data = t) %>%
        predict(newdata = data.frame(x = max(t$x)))
      , max = max(t$x)
    )}) %>%
  bind_rows


CesAplot <- ggplot(CesAexpr, aes(x=DPA,y=expression,group=homoeolog)) +
	geom_smooth(aes(color = stage),
            alpha = 0.2) +
	geom_dl(aes(label = CesA), 
		method=list("last.points", cex=0.7, hjust=-0.1))  +
	xlim(5,30) +
	ylim(0,15) +
 	theme_bw() +
	theme(axis.text.x = element_text(angle = 90))+
	facet_wrap(~stage, ncol=1)

ggsave("CesAExpression.jpg", plot=CesAplot, width = 10, height = 15, units = "in")




CesAplotH <- ggplot(CesAexpr, aes(x=DPA,y=expression,group=homoeolog)) +
	geom_smooth(aes(color = stage),
            alpha = 0.2) +
	xlim(5,30) +
	ylim(0,15) +
 	theme_bw() +
	theme(axis.text.x = element_text(angle = 90))+
	facet_wrap(~homoeolog, ncol=5)

ggsave("CesAExpressionHomoeolog.jpg", plot=CesAplotH, width = 10, height = 10, units = "in")


