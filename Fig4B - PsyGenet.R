library("psygenet2r")

###################

setwd('/path/to/my/directory/')

allresults <- as.data.frame(read.table(file=file.path(currdir, 'All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt'), header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

baseMean <- lapply(names(allresults), function(ch) grep('baseMean', ch))
baseMean <- allresults[, which(baseMean==1)] 
minimums <- unlist(apply(baseMean, 1, min))  

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

cats <- allresults[, c(241:288)]

################### 
# PsyGenet Database, in text form
total <- read.table(file='psygenet_v02.txt', sep='\t', header=TRUE, stringsAsFactors=F)
PsyDisorder <- unique(total$PsychiatricDisorder)
Diseasename <- total$DiseaseName

# Incorporating a human-rat gene mapping database provided by SynGo.
mappedgenes <- read.csv(file='psygenet_v02_idmap.csv', stringsAsFactors=F)
mapped <- c()
for(i in 1:nrow(total)){
	mapped <- c(mapped, mappedgenes$symbol[match(total$Gene_Symbol[i], mappedgenes$query)])
}
total <- cbind(total, mapped)

###################
# Specify down-regulated, up-regulated or all DEGs
classification = c('All', 'Down', 'Up')

for(class in classification){
	masteralldiseases <- list()	
	
	for(h in 1:length(PsyDisorder)){ 
		listofgenesall <- lapply(total$PsychiatricDisorder, function(ch) grep(PsyDisorder[h], ch))
		listofgenesall <- total$mapped[which(listofgenesall==1)]
		listofgenesall <- unique(listofgenesall)
		
		listofgenes <- intersect(listofgenesall, row.names(allresults))  
		
		master <- matrix(0, ncol(cats), 10)
		
		for(i in 1:ncol(cats)){
			frame <- allresults[, c(((i-1)*5+3),((i-1)*5+4), 240+i)]
			if(class=='All'){
				frame2 <- frame[frame[[3]]!=0, ]
			} else if(class=='Down'){
				frame2 <- frame[frame[[3]]<0, ]
			} else if(class=='Up'){
				frame2 <- frame[frame[[3]]>0, ]
			}
			  		
			genenames <- row.names(frame2)  #### After mapping with SynGO, all gene names are concordant.
			
			found <- intersect(genenames, listofgenes)
			
			fisherstable <- matrix(0, 2, 2)
			fisherstable[1,1] <- length(found)
			fisherstable[1,2] <- length(listofgenes) - length(found)
			fisherstable[2,1] <- length(genenames) - length(found)
			fisherstable[2,2] <- length(row.names(frame)) - length(genenames) - fisherstable[1,2]
			
			fisherstest <- fisher.test(fisherstable)
			
			master[i, 1] <- -log(fisherstest$p.value, 10)
			master[i, 2:6] <- as.numeric(unlist(fisherstest)[1:5])
			master[i, 7:8] <- fisherstable[1, ]
			master[i, 9:10] <- fisherstable[2, ]		
		}
		
		master <- as.data.frame(master)
		row.names(master) <- names(cats)
		names(master) <- c('(-log10 p-value)', attr(unlist(fisherstest), "names")[1:5], 'Intersect', 'In Pathway', 'DE-Intersect', 'Rest')
		
		name <- PsyDisorder[h]; name <- gsub('/', '_', name) # simplify PsyDisorder name
		masteralldiseases[[name]] <- master
	}
	
	######################
	# PERFORM BH CORRECTION
	for(disorder in 1:length(masteralldiseases)){
		data <- masteralldiseases[[disorder]]
		datapvalues <- data$p.value
		datapvaluesBH <- p.adjust(datapvalues, method='BH')
		data[['adjusted -log10 p-value']] <- -log10(datapvaluesBH)
		data[['adjusted p.value']] <- datapvaluesBH
		data <- data[, c(11,12,1:10)]
		masteralldiseases[[disorder]] <- data	
		write.csv(as.data.frame(data), file=paste0('Psygenet_', class, '_', names(masteralldiseases)[disorder], '_Fishers analysis-BH-corr.csv'))
	}

}