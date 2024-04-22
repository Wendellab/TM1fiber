library(WGCNA)
library(RColorBrewer)
library(data.table)
library(DESeq2)
library(tidyverse)

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



#####################
## 2. Choose Power ##
#####################

load("04.R-01-dataInput.filtered.RData")


### Choose a set of soft-thresholding powers ### 

# powers to test
powers = c(c(1:10), seq(from = 12, to=40, by=2))


# Call the network topology analysis function 
powerTable = pickSoftThreshold(datExpr, powerVector=powers, verbose = 2, networkType = "signed")[[2]]      
collectGarbage()
    



### Plot the results ###
colors=brewer.pal(5,"Set1")

# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)

for (col in 1:length(plotCols)) {
        ylim[1, col] = min(ylim[1, col], powerTable[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTable[, plotCols[col]], na.rm = TRUE);
}

    
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)

pdf("out-05.s2.ChooseSoftThresholdPower.filtered.pdf" )
    par(mfcol = c(2,2));
    par(mar = c(4.2, 4.2 , 2.2, 0.5))
    cex1 = 0.7;


  for (col in 1:length(plotCols)) {
   plot(powerTable[,1], -sign(powerTable[,3])*powerTable[,2], xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col]);
            addGrid()
   
    if (col==1)
    {
        text(powerTable[,1], -sign(powerTable[,3])*powerTable[,2],
        	labels=powers,cex=cex1,col=colors[2]);
    } else
    	text(powerTable[,1], powerTable[,plotCols[col]], labels=powers,cex=cex1,col=colors[2]);
    }

  
dev.off()

save(powerTable, file = "05.R-02-choosePower.filtered.RData")

# Inspect "s2.ChooseSoftThresholdPower_signed.pdf".
# Basically, the default power is 6 for unsigned network, 12 for signed network; 
# if the observed power choice is smaller than default, I will choose the observed power, 
# otherwise use default
### 10 is good here
