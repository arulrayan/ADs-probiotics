library('RRHO2')

#####################
setwd('/path/to/my/directory/')

allresults <- as.data.frame(read.table(file=file.path(currdir, 'All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt'), header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

padj <- lapply(names(allresults), function(ch) grep('padj', ch))
padj <- allresults[, which(padj==1)]

log2FCF <- log2FC
padjF <- padj  ## NECESSARY BECAUSE PADJ IS CALLED BELOW

cats <- allresults[, c(241:288)]
catsF <- cats
allsums <- allresults[, 289]
allabssums <- allresults[, 290]

#####################

setwd('/path/to/my/directory/Girgenti')
files <- list.files(getwd(), pattern='full.csv')

#####################
master <- NULL
for(file in files){
	if(grepl('dacc', file)){
		regions <- c(34:38)
	}
	if(grepl('sgpfc', file)){
		regions <- c(15:19)
	}
	data <- read.csv(file, stringsAsFactors=F)
	genes <- data$Genename
	logFC <- data[[3]]
	padj <- data[[4]] #### THE VALUE PROVIDED IS ACTUALLY JUST THE PADJ
	padjtransform <- -log10(padj)*logFC/abs(logFC)
	
	available <- match(toupper(genes), toupper(rownames(catsF)))
	availableF <- available[!is.na(available)]
	availableF <- unique(availableF)
	availablegenes <- rownames(catsF)[availableF]
	matched <- match(toupper(availablegenes), toupper(genes))
	matchedF <- matched[!is.na(matched)]
	matchedF <- unique(matchedF)
	
	for(region in regions){
		regiondata <- log2FCF[availableF, region]
		logFCdata <- logFC[matchedF]
		padjdata <- padj[matchedF]
		
		padjtransformF <- padjtransform[matchedF]
		padjtransformF[is.na(padjtransformF)] <- 0
		
		padjF2 <- padjF[availableF, region]
		padjF2transform <- -log10(padjF2)*regiondata/abs(regiondata)
		padjF2transform[is.na(padjF2transform)] <- 0
		
		if(length(padjF2transform)>=0){
			Ext <- data.frame(Genes=availablegenes, DDE=padjtransformF, stringsAsFactors=F)
			Int <- data.frame(Genes=availablegenes, DDE=padjF2transform, stringsAsFactors=F)
			RRHO_obj <-  RRHO2_initialize(Ext, Int, labels = c(file, names(catsF)[region]), log10.ind=TRUE, multipleTesting='BH', boundary=0.02)
			
			filename <- paste("RRHOMap_padj_", names(catsF)[region], '_', substring(file, 1, nchar(file)-4), ".tiff", sep = "")
        	tiff(filename = filename, width = 6, height = 6, units = "in", res = 150)
			RRHO2_heatmap(RRHO_obj) 
			dev.off()			
		}
	}
}
	
