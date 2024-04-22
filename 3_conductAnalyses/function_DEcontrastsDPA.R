###########################################
##### partitioned homoeolog contrasts #####
###########################################

contrastTable <- data.frame(contrast=character(),
	DE=integer(),DEup=integer(),DEdown=integer(),
	DEA=integer(),DEAup=integer(),DEAdown=integer(),
	DED=integer(),DEDup=integer(),DEDdown=integer(), 
	stringsAsFactors=F)
	
SigFoldChange <- data.frame(target_id=rownames(dds))

DPAcontrast <- function (DPA2, DPA1) {
	DPAres <- results(dds,contrast=c("DPA",paste0(DPA2,"DPA"),paste0(DPA1,"DPA")))
	name <- paste0("contrast",DPA2,"v",DPA1)	

	DE       <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% nrow()
	DEup     <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEdown   <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()
	
	DEA      <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% nrow()
	DEAup    <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEAdown  <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()

	DED      <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% nrow()
	DEDup    <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEDdown  <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()

    contrastTable[nrow(contrastTable)+1,] <- rbind(paste0(DPA2," v ", DPA1),DE,DEup,DEdown,DEA,DEAup,DEAdown,DED,DEDup,DEDdown)        
	assign("contrastTable",contrastTable,envir = globalenv())

	newSFC <- DPAres %>%
		as.data.frame %>%
		filter(padj<0.05) %>%
		select("log2FoldChange") %>%
		rownames_to_column("target_id") %>% 
		dplyr::rename(!!name := log2FoldChange) %>%
		left_join(SigFoldChange, .)

	assign("SigFoldChange",newSFC,envir = globalenv())	

	write.table(DPAres,file=paste0("DEcontrasts/out-02.DEGresults.",DPA2,"v",DPA1,".tsv"),quote=F, sep="\t")
}




######################################
##### summed homoeolog contrasts #####
######################################

ScontrastTable <- data.frame(contrast=character(),
	DE=integer(),DEup=integer(),DEdown=integer(), 
	stringsAsFactors=F)

SDPAcontrast <- function (DPA2, DPA1) {
	DPAres <- results(Sdds,contrast=c("DPA",paste0(DPA2,"DPA"),paste0(DPA1,"DPA")))
	
	DE       <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% nrow()
	DEup     <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEdown   <- as.data.frame(DPAres) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()
	
	DEA      <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% nrow()
	DEAup    <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEAdown  <- as.data.frame(DPAres) %>% filter(grepl("\\.A", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()

	DED      <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% nrow()
	DEDup    <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% nrow()
	DEDdown  <- as.data.frame(DPAres) %>% filter(grepl("\\.D", rownames(DPAres))) %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% nrow()

    ScontrastTable[nrow(ScontrastTable)+1,] <- rbind(paste0(DPA2," v ", DPA1),DE,DEup,DEdown)        
	assign("ScontrastTable",ScontrastTable,envir = globalenv())
	
	write.table(DPAres,file=paste0("DEcontrasts/out-02.DEGresults.",DPA2,"v",DPA1,".summed.tsv"),quote=F, sep="\t")
}