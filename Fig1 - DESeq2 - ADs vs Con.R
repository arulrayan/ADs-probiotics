library("DESeq2")
library('dplyr')

setwd("/path/to/my/directory/Bulk RNAseq/Trmts_vs_Con")

filelistraw <- list.files(getwd(), pattern = '_.txt')
filenumbers <- c(1:length(filelistraw))
filelistrawlib <- list.files(getwd(), pattern = '_Lib.txt')

print(filelistrawlib)

for(h in 1:length(filelistraw)){
	dataraw <- as.data.frame(read.table(filelistraw[h], header=TRUE, sep='\t', stringsAsFactors=F, row.names='gene_names'))	
	datarawname <- substring(filelistraw[h], unlist(gregexpr('sortedby_', filelistraw[h]))+9, unlist(gregexpr('_.txt', filelistraw[h]))-1)

	rawlibindex <- unlist(gregexpr(datarawname, filelistrawlib)); rawlibindex=which(rawlibindex>0)
	datalibraw <- as.data.frame(read.table(filelistrawlib[rawlibindex], header=TRUE, sep='\t'))
	datalibraw$libraryID <- as.character(datalibraw$libraryID)

	datalibraw$drug <- as.factor(datalibraw$drug)
	drugsraw <- levels(datalibraw$drug)
	
	ddsbatch <- DESeqDataSetFromMatrix(countData = dataraw, colData = datalibraw, design = ~drug + Replicate)
	print(ddsbatch)
	keepbatch <- rowSums(counts(ddsbatch)) >= 5*(ncol(dataraw))
	ddsbatch <- ddsbatch[keepbatch,]
	ddsbatch$drug <- relevel(ddsbatch$drug, ref = 'Control')
	ddsDESeqbatch <- DESeq(ddsbatch)

	masternormcounts <- NULL
	normdds <- estimateSizeFactors(ddsDESeqbatch)
	normcounts <- counts(normdds, normalized=TRUE)
	masternormcounts <- as.data.frame(normcounts)
	
	masterresultsbatch <- list()
	masterresultsmatrixbatch <- matrix(0, nrow=length(drugsraw), ncol=length(drugsraw))
	
	masterresbatch <- as.data.frame(mcols(ddsDESeqbatch))
	masterresbatch <- masterresbatch[,c(1,2,4,7)]
	
	masterres <- results(ddsDESeqbatch)
	
	#FOR EACH RESULTSNAMES, ADD 1)LOG2FOLDCHANGE, 2)LFCSE, 3)PVALUE AND 4)PADJ, IN THAT ORDER.	
	for(j in 2:(length(drugsraw)-1)){  ## DON'T USE THE PROBIOTICS COMPARISON
		 #orders, by p-value, the results for that particular comparison and saves.
		 resbatchraw <- results(ddsDESeqbatch, name=resultsNames(ddsDESeqbatch)[j])
		 resbatchLFC <- lfcShrink(ddsDESeqbatch, coef=resultsNames(ddsDESeqbatch)[j], type='apeglm')
		 resbatch <- as.data.frame(resbatchLFC)
		 resbatch <- resbatch[, c(2,3,4,5)]
		 drugnamebatch <- substring(resultsNames(ddsDESeqbatch)[j], 6)
		 locofunderscorebatch <- gregexpr("_", drugnamebatch)
		 if (locofunderscorebatch[[1]][1]>3){
		 	drugnamebatch <- substring(drugnamebatch, 1, 3)
		 }	else {
		 	drugnamebatch <- substring(drugnamebatch, 1, 2)
		 }		 
		 names(resbatch) <- c( paste(drugnamebatch, '_log2FC', sep=""), paste(drugnamebatch, '_lFCSE', sep=""), paste(drugnamebatch, '_pvalue', sep=""), paste(drugnamebatch, '_padj', sep="") )
		 masterresbatch <- cbind(masterresbatch, resbatch)
	}
	write.csv(as.data.frame(masterresbatch), file=paste(substring(filelistraw[h], 1, nchar(filelistraw[h])-5),"-DESeq2w-R-corr-APEGLM-ILCsLC-MR.csv",sep=""))
	write.table(as.data.frame(masternormcounts), file=paste(substring(filelistraw[h], 1, nchar(filelistraw[h])-5),"-DESeq2w-R-corr-APEGLM-masternormcounts.txt",sep=""), quote=F, sep='\t')

}
	 
