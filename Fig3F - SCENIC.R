library(pheatmap)
library(SCENIC)
library(WGCNA)
library("dplyr")
library("tibble")

###################
# function to transform log2FC values
log2FCtransform <- function(listofvalues){   ##listofvalues must be a list
	elistofvalues <- unlist(lapply(listofvalues, exp))
	return(elistofvalues/(elistofvalues+1))
}
###################
# SCENIC databases
dbDir="/Library/Frameworks/R.framework/Versions/3.6/Resources/library/SCENIC"
dbs='mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather'

###################
currdir <- '/path/to/my/directory/'

allresults <- as.data.frame(read.table(file=file.path(currdir, 'All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt'), header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

baseMean <- lapply(names(allresults), function(ch) grep('baseMean', ch))
baseMean <- allresults[, which(baseMean==1)] 
minimums <- unlist(apply(baseMean, 1, min))  

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

###################
samplenames <- substring(names(log2FC), 1, nchar(names(log2FC))-7)
masterlog2FCtransform <- c()
for(i in 1:ncol(log2FC)){
	try <- log2FCtransform(log2FC[[i]])
	masterlog2FCtransform <- cbind(masterlog2FCtransform, try)
}
attr(masterlog2FCtransform, "dimnames")[[1]] <- rownames(log2FC)
attr(masterlog2FCtransform, "dimnames")[[2]] <- samplenames

###################
outdir <- dir.create(file.path(currdir, 'SCENIC'), showWarnings=F)
setwd(file.path(currdir, 'SCENIC'))

avgexprMatF <- masterlog2FCtransform
avgscenicOptionsF <- initializeScenic(org="mgi", dbDir=dbDir, dbs=dbs)
saveRDS(avgscenicOptionsF, file="int/scenicOptions.Rds")

runCorrelation(avgexprMatF, avgscenicOptionsF)
runGenie3(avgexprMatF, avgscenicOptionsF)

avgscenicOptionsF2 <- avgscenicOptionsF
saveRDS(avgscenicOptionsF2, file='int/scenicOptionsbackup.Rds')

avgscenicOptionsF2@settings$nCores <- 4
avgscenicOptionsF2@settings$seed <- 123
runSCENIC_1_coexNetwork2modules(avgscenicOptionsF2)
runSCENIC_2_createRegulons(avgscenicOptionsF2)
runSCENIC_3_scoreCells(avgscenicOptionsF2, avgexprMatF)

################

avgregulonAUC <- readRDS('int/3.4_regulonAUC.Rds')
avgregulonAUC <- avgregulonAUC[onlyNonDuplicatedExtended(rownames(avgregulonAUC)),]
avgdata <- avgregulonAUC@assays@data@listData$AUC

avgdataF <- as.data.frame(avgdata)
write.table(avgdataF, 'Master Results, all nonDupExtend regulon AUCs.txt', sep='\t', quote=F)