library("DESeq2")
library('dplyr')

setwd("/path/to/my/directory/Bulk RNAseq/Pro_vs_Mal")

datafiles <- list.files(getwd(), pattern = '_.txt')
libfiles <- list.files(getwd(), pattern = '_Lib.txt')

for(i in 1:length(datafiles)){
	filename <- datafiles[i]
	cts <- as.matrix(read.csv(datafiles[i], sep="\t", row.names="gene_names"))
	coldata <- read.csv(libfiles[i], sep="\t", row.names=1)
	
	coldatadrug <- coldata[coldata$drug!='Maltodextrin',]
	drugunique <- unique(coldatadrug$region)
	coldatacontrol <- coldata[coldata$drug=='Maltodextrin',]
 	controlunique <- unique(coldatacontrol$region)
	regions <- intersect(drugunique, controlunique)
	
	##################
	##################
	
	##INSERT CODE THAT: SEARCHES COLDATA FOR BRAIN REGION (USE FILTER), THEN RETURNS SAMPLE IDS, AND KEEPS SAMPLE IDS IN CTS. THEN DO DESEQ.
	masterres <- NULL
	masternormcounts <- NULL
	for (brainregion in 1:length(regions)){
		coldataTF <- coldata$region==regions[brainregion]
		coldatasubset <- coldata[coldataTF,]
		ctssubset <- cts[, coldataTF]
		
		### BY REGION:: #FILTER BY AVERAGE MEAN CTS (at least 5)
		sumcts <- apply(ctssubset, 1, sum)
		ctssubset <- ctssubset[sumcts>=(5*dim(ctssubset)[2]), ]
		
		dds <- DESeqDataSetFromMatrix(countData = ctssubset, colData = coldatasubset, design = ~drug + Replicate)
	    print(dds)
	    numcontrols <- sum(coldatasubset$drug=='Maltodextrin')
	    totsamples <- length(names(coldatasubset))
		 
		dds$drug <- relevel(dds$drug, ref = "Maltodextrin")
		ddsDESeq <- DESeq(dds)
		
		resLFC <- lfcShrink(ddsDESeq, coef=resultsNames(ddsDESeq)[2], type='apeglm')
		
		res <- as.data.frame(resLFC)
		drugname <- substring(resultsNames(ddsDESeq)[2], 6)
		locofunderscore <- gregexpr("_", drugname)
		if (locofunderscore[[1]][1]>3){
			drugname <- substring(drugname, 1, 3)
		}	else {
			drugname <- substring(drugname, 1, 2)
		}		
		drugregion <- paste(drugname, regions[brainregion], sep='_') 
		names(res) <- c(paste(drugregion, '_baseMean', sep=''), paste(drugregion, '_log2FC', sep=""), paste(drugregion, '_lFCSE', sep=""), paste(drugregion, '_pvalue', sep=""), paste(drugregion, '_padj', sep="") )
							
		normdds <- estimateSizeFactors(ddsDESeq)
		normcounts <- as.data.frame(counts(normdds, normalized=TRUE))
		
		if(brainregion==1){
			masterres <- res
			masternormcounts <- as.data.frame(normcounts)
		} else {
			masterres <- merge(masterres, res, by=0, all=TRUE)
			rownames(masterres) <- masterres$Row.names; 
			masterres$Row.names <- NULL
			
			masternormcounts <- merge(masternormcounts, normcounts, by=0, all=TRUE)	
			rownames(masternormcounts) <- masternormcounts$Row.names; 
			masternormcounts$Row.names <- NULL	
		}		
	}
	
	write.csv(as.data.frame(masterres), file=paste(substring(datafiles[i], 1, nchar(datafiles[i])-5),"-DESeq2w-R-corr-APEGLM-masterres.csv",sep=""))
	write.csv(as.data.frame(masternormcounts), file=paste(substring(datafiles[i], 1, nchar(datafiles[i])-5),"-DESeq2w-R-corr-APEGLM-masternormcounts.csv",sep=""))
	
}