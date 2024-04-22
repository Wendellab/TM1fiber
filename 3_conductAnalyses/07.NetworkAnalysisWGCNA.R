library(WGCNA)
library(RColorBrewer)
library(scatterplot3d)
library(data.table)
library(DESeq2)
library(tidyverse)
library(flashClust)
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


#########################
## 4. Network Analysis ##
#########################

source('multiplot.r', chdir = TRUE)
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('WGCNA_missingFun.r', chdir = TRUE)

load(file = "R-01-dataInput.filtered.RData") 
load(file = "04.R-01-dataInput.filtered.RData")
load(file = "06.R-03-buildNetwork.RData")


### write gene with corresponding module assignment and annotation ###
aa <- fread(file="aa.annotation.data.tsv", header=T)

me <- data.frame(homoeolog = colnames(datExpr), ME_rld = net$colors)

me$gene <- gsub("[.]1[.][A|D]$","", me$homoeolog)
me$sub <- gsub(".*[.]","", me$homoeolog)
table(me$ME_rld,me$sub)
dim(me <- merge(me, aa, all.x=TRUE, by="gene"))
fwrite(me, file="out-07.s3.module_annotation.txt", row.names=FALSE,sep="\t")

MEs <- net$MEs
Nmodules <- dim(net$MEs)[2]
  
print(net$TOMFiles)
print(paste("Number of modules in network is ", Nmodules, sep=""))

pval <- as.data.frame((apply(MEs,2,function(x){round(anova(aov(x~coldata$DPA) )$"Pr(>F)"[1],4)})))
names(pval) <-"sample"
pvals <- ifelse(pval<0.05,"*","")



### add number of genes per module ###
nn <- as.data.frame(table(net$colors))
names(nn) <-c ("module","geneN")
rownames(nn) <- paste0("ME",nn$module)
pvals <- pvals[rownames(nn),]
nsig <- cbind(nn,pvals)
write.table(nsig, file="out-07.s4.filtered.moduleSignificance.txt", row.names=FALSE,sep="\t")
    


### make multicolor barplots ###
pdf("out-07.s5.filtered.modules.pdf")
plots <- list()  # new empty list

for(me in 0:(Nmodules-1)) {
   which.module=paste("ME",me,sep="")
   module.color=labels2colors(me)

   # heatmap
   	par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
	   plotMat(t(scale(datExpr[,net$colors==me ]) ),
           nrgcols=30,rlabels=T,rcols=module.color,
           main=paste0(which.module, ": ", module.color, " ",nsig[which.module,"*"]), cex.main=2)

   # barplot
   	par(mar=c(5, 4.2, 0, 0.7))
	   barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
	   ylab="eigengene expression",xlab="sample", names.arg=coldata$DPA,las=2,cex.names=0.6)

   # MEplots, anova
   	df <- data.frame(ME=MEs[,which.module], module = which.module,dpa=coldata$DPA,sample=coldata$DPA )
	   fit <- aov(ME ~sample, df)
	   dfc <- summarySE(df, measurevar="ME", groupvars=c("sample","module","dpa"))

	plots[[me+1]] <- ggplot(dfc, aes(x=sample, y=ME, fill = dpa)) + 
		geom_bar(stat="identity",position=position_dodge(), color="black", size=0.3) + 
			geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.3,position=position_dodge(0.9)) + 
			ggtitle(paste0(which.module," ",module.color),
				subtitle= paste0("genes: ", nsig[which.module,"geneN"], ", P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") ) + 
			theme_bw() + 
			theme(plot.title = element_text(size=11, hjust = 0.5), 
				plot.subtitle = element_text(size=9, hjust = 0.5), 
				axis.text.x = element_text(size = 6, angle = 65, hjust=1), 
				legend.position = "none")
}

for(page in 1:ceiling(Nmodules/9)) {
    if(Nmodules>(9*page)) {  
	multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
     else {  
	multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
}

plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

dev.off()
# end the big plot


pdf("out-07.s5.filtered.modules.multiplot.pdf")

for(page in 1:ceiling(Nmodules/9)) {
    if(Nmodules>(9*page)) {  
	multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
     else {  
	multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
}

dev.off()


save(Nmodules, MEs, file = "07.R-04-networkTopology.RData")

######################
### module heatmap ###
######################

  data <- cbind(MEs, coldata) %>%
    gather(starts_with("ME"), key="module", value="eigengene") %>%
    group_by(module, species, DPA) %>%
    summarise(SE = sd(eigengene)/sqrt(n()), eigengene=mean(eigengene)) %>%
    ungroup() %>%
    mutate(module=factor(module, rev(sprintf("ME%d", 0:17)))) %>%
    mutate(SE = ifelse(SE > 0.1, 0.1, SE)) %>%
    mutate(DPA = as.numeric(substring(DPA,1,2))) 

  module.info <- nsig %>%
    select(-'module') %>%
    rownames_to_column('module')

    
  heat.plot <- ggplot(data, aes(DPA, module, color=eigengene, size=SE)) +
    geom_point(shape=15) +
    geom_text(aes(x=25, y=module, label=scales::comma(geneN)), module.info, size=4, color='black') +
    scale_color_distiller(type='div', palette='RdBu', direction=-1, limits=c(-0.35, 0.35)) +
    scale_size_continuous(trans='reverse', range=c(3, 10)) +
    scale_x_continuous(breaks=c(6,10,15,20,24)) +
    scale_y_discrete(labels=rev(sprintf("%s%s",
                                        module.info$pvals,
                                        module.info$module))) +
    coord_equal() +
  theme_minimal() +
    theme(axis.title.y=element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
  #        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.text.y=element_text(face=ifelse(rev(module.info$pval)=='*',
                                               'bold', 'plain'),
                                   size = 10)
          ) 

  heat.plot +
    annotate(geom='segment', x=9.5, xend=9.5, y=0.5, yend=18.5)+
    annotate(geom='segment', x=16.5, xend=16.5, y=0.5, yend=18.5)+
    annotate(geom='errorbar', y=19, xmin=8.5, xmax=12.5, width=0.5) +
    annotate(geom='label', label="Rapid Elongation",  y=19, x=10.5) +
    annotate(geom='errorbar', y=19, xmin=15.5, xmax=20.5, width=0.5) +
    annotate(geom='label', label="Transition",  y=19, x=18) +
    annotate(geom='errorbar', y=20, xmin=13.5, xmax=24.5, width=0.5) +
    annotate(geom='label', label="Secondary Wall Synthesis",  y=20, x=18.5) +
    annotate(geom='errorbar', y=21, xmin=5.5, xmax=20.5, width=0.5) +
    annotate(geom='label', label="Primary Wall Synthesis",  y=21, x=13) +
    annotate(geom='errorbar', color = "white", y=21.5, xmin=5.5, xmax=20.5, width=0.5)


  ggsave("out-07.ModuleHeatmap.png", bg='white', dpi=600)


