library(BRETIGEA, quietly=T)

###################

currdir <- '/path/to/my/directory/'

allresults <- as.data.frame(read.table(file=file.path(currdir, 'All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt'), header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

baseMean <- lapply(names(allresults), function(ch) grep('baseMean', ch))
baseMean <- allresults[, which(baseMean==1)] 
minimums <- unlist(apply(baseMean, 1, min))  

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 
cats <- allresults[, c(241:288)]

###################
# Load BRETIGEA database

mousebrainmarkers <- markers_df_mouse_brain
celltypes <- unique(mousebrainmarkers$cell)

###################
# Specify down-regulated, up-regulated or all DEGs
classification = c('All', 'Down', 'Up')

for(class in classification){
	masterallcelltypes <- list()
	
	for(celltype in celltypes){  # for six cell types in the database
		
		master <- matrix(0, ncol(cats), 10)
		
		for(sample in 1:ncol(cats)){  # for 48 RTPs
			totalgenesinsampleBM <- allresults[, (sample-1)*5+1]
			totalgenesinsamplenames <- rownames(allresults)[totalgenesinsampleBM>0]
			totalgenesinsample <- sum(totalgenesinsampleBM>0)
			
			genesavailable <- intersect(totalgenesinsamplenames, mousebrainmarkers$markers[mousebrainmarkers$cell==celltype])
			
			if(class=='All'){
				genes <- rownames(allresults)[cats[sample]!=0] ##by default, these are already categorized, so they'll have no problems fitting the criteria.
			} else if(class=='Down'){
				genes <- rownames(allresults)[cats[sample]<0]
			} else if(class=='Up'){
				genes <- rownames(allresults)[cats[sample]>0]
			}
			found <- intersect(genes, genesavailable)   
			
			fisherstable <- matrix(0, 2, 2)
			fisherstable[1,1] <- length(found)
			fisherstable[1,2] <- length(genesavailable) - length(found)
			fisherstable[2,1] <- length(genes) - length(found)
			fisherstable[2,2] <- totalgenesinsample - length(genes) - fisherstable[1,2]
			
			fisherstest <- fisher.test(fisherstable)  # construct Fisher's exact test
			
			master[sample, 1] <- -log(fisherstest$p.value, 10)
			master[sample, 2:6] <- as.numeric(unlist(fisherstest)[1:5])
			master[sample, 7:8] <- fisherstable[1, ]
			master[sample, 9:10] <- fisherstable[2, ]
		}
		
		master <- as.data.frame(master)
		row.names(master) <- names(cats)
		names(master) <- c('(-log10 p-value)', attr(unlist(fisherstest), "names")[1:5], 'Intersect', 'In Pathway', 'DE-Intersect', 'Rest')
		masterallcelltypes[[celltype]] <- master	
	}
	
	######################
	# PERFORM BH CORRECTION
	for(celltype in 1:length(masterallcelltypes)){		
		data <- masterallcelltypes[[celltype]]
		datapvalues <- data$p.value
		datapvaluesBH <- p.adjust(datapvalues, method='BH')
		data[['adjusted -log10 p-value']] <- -log10(datapvaluesBH)
		data[['adjusted p.value']] <- datapvaluesBH
		data <- data[, c(11,12,1:10)]
		masterallcelltypes[[celltype]] <- data
		write.csv(as.data.frame(masterallcelltypes[[celltype]]), file=paste0('BRETIGEA_', class, '_', names(masterallcelltypes)[celltype], '_Fishers analysis-BH-corr.csv'))
	}

}	
