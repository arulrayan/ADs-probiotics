library(Seurat)
library(DESeq2)

setwd('/path/to/my/directory/')

ILCLMPFS <- readRDS('ILCL.rds')
ILCRMPFS <- readRDS('ILCR.rds')
PLCLMPFS <- readRDS('PLCL.rds')
PLCRMPFS <- readRDS('PLCR.rds')

################
################

param <- 'Banksy_lam0.15_k20_res0.5'

#### Count the number of spots per cluster, across all 24 pseudobulk samples for ILC and PLC separately. Only use clusters with >=3 spots, on average, across the 24 pseudobulk samples. This yields clusters 3,4,6,7,8. 3 and 4 are merged due to transcriptomic similarity.
ILCLmeta <- ILCLMPFS@meta.data
ILCRmeta <- ILCRMPFS@meta.data
ILCclusts <- unique(c(ILCLmeta[[param]], ILCRmeta[[param]]))
ILCclusts <- sort(ILCclusts)
ILCmat <- matrix(0, nrow=length(unique(ILCLmeta$SampleType))*2, ncol=length(unique(c(ILCLmeta[[param]], ILCRmeta[[param]]))))
for(i in 1:(nrow(ILCmat)/2)){
	dataL <- ILCLmeta[ILCLmeta$SampleType==unique(ILCLmeta$SampleType)[i], ]
	dataR <- ILCRmeta[ILCRmeta$SampleType==unique(ILCRmeta$SampleType)[i], ]
	for(j in 1:ncol(ILCmat)){
		ILCmat[(2*i-1),j] <- sum(dataL[[param]]==ILCclusts[j])
		ILCmat[2*i,j] <- sum(dataR[[param]]==ILCclusts[j])		
	}
}
ILCmat <- as.data.frame(ILCmat, row.names=sort(c(paste0(unique(ILCLmeta$SampleType),'_L'), paste0(unique(ILCRmeta$SampleType),'_R')))) 
names(ILCmat) <- paste0('Clust_', ILCclusts)
ILCcluststouse <- ILCclusts[unlist(apply(ILCmat,2,sum))>=3*nrow(ILCmat)]

PLCLmeta <- PLCLMPFS@meta.data
PLCRmeta <- PLCRMPFS@meta.data
PLCclusts <- unique(c(PLCLmeta[[param]], PLCRmeta[[param]]))
PLCclusts <- sort(PLCclusts)
PLCmat <- matrix(0, nrow=length(unique(PLCLmeta$SampleType))*2, ncol=length(unique(c(PLCLmeta[[param]], PLCRmeta[[param]]))))
for(i in 1:(nrow(PLCmat)/2)){
	dataL <- PLCLmeta[PLCLmeta$SampleType==unique(PLCLmeta$SampleType)[i], ]
	dataR <- PLCRmeta[PLCRmeta$SampleType==unique(PLCRmeta$SampleType)[i], ]
	for(j in 1:ncol(PLCmat)){
		PLCmat[(2*i-1),j] <- sum(dataL[[param]]==PLCclusts[j])
		PLCmat[2*i,j] <- sum(dataR[[param]]==PLCclusts[j])		
	}
}
PLCmat <- as.data.frame(PLCmat, row.names=sort(c(paste0(unique(PLCLmeta$SampleType),'_L'), paste0(unique(PLCRmeta$SampleType),'_R')))) 
names(PLCmat) <- paste0('Clust_', PLCclusts)
PLCcluststouse <- PLCclusts[unlist(apply(PLCmat,2,sum))>=3*nrow(PLCmat)]

################
################
#### For the clusters, create pseudobulks. If cluster is 3, omit; if cluster is 4, count both cluster 4 and cluster 3 spots - this merges clusters 3 & 4.
ILCMPFSLRrawPB <- list()
for(cluster in ILCcluststouse){
	print(cluster)
	for(sample in unique(ILCLmeta$SampleType)){
		subsetted <- subset(ILCLMPFS, subset=SampleType==sample)
		rawdata <- subsetted$Spatial@data
		rawdatamat <- as.matrix(rawdata)
		if(cluster==3){ 
			cells <- NULL
		} else if(cluster==4){   ## THIS MERGES CLUSTERS 3 & 4
			cells <- subsetted@meta.data[[param]]==cluster | subsetted@meta.data[[param]]==cluster-1
		} else {
			cells <- subsetted@meta.data[[param]]==cluster
		}
		if(!is.null(cells)){
			if(sum(cells)>1){
				rawdatamatF <- rawdatamat[, cells]
				if(sum(cells)>1){ rawPB <- rowSums(rawdatamatF)} 
				else { rawPB <- rawdatamatF }
				ILCMPFSLRrawPB[[paste0(sample,'_L_clust_', cluster)]] <- rawPB
			}
		}
	}
	for(sample in unique(ILCRmeta$SampleType)){
		subsetted <- subset(ILCRMPFS, subset=SampleType==sample)
		rawdata <- subsetted$Spatial@data
		rawdatamat <- as.matrix(rawdata)
		if(cluster==3){ 
			cells <- NULL
		} else if(cluster==4){    ## THIS MERGES CLUSTERS 3 & 4
			cells <- subsetted@meta.data[[param]]==cluster | subsetted@meta.data[[param]]==cluster-1
		} else {
			cells <- subsetted@meta.data[[param]]==cluster
		}
		if(!is.null(cells)){
			if(sum(cells)>1){
				rawdatamatF <- rawdatamat[, cells]
				if(sum(cells)>1){ rawPB <- rowSums(rawdatamatF)} 
				else { rawPB <- rawdatamatF }
				ILCMPFSLRrawPB[[paste0(sample,'_R_clust_', cluster)]] <- rawPB
			}
		}
	}
}

#### PLC
PLCMPFSLRrawPB <- list()

for(cluster in PLCcluststouse){
	print(cluster)
	for(sample in unique(PLCLmeta$SampleType)){
		subsetted <- subset(PLCLMPFS, subset=SampleType==sample)
		rawdata <- subsetted$Spatial@data
		rawdatamat <- as.matrix(rawdata)
		if(cluster==3){ 
			cells <- NULL
		} else if(cluster==4){
			cells <- subsetted@meta.data[[param]]==cluster | subsetted@meta.data[[param]]==cluster-1
		} else {
			cells <- subsetted@meta.data[[param]]==cluster
		}
		if(!is.null(cells)){
			if(sum(cells)>1){
				rawdatamatF <- rawdatamat[, cells]
				if(sum(cells)>1){ rawPB <- rowSums(rawdatamatF)} 
				else { rawPB <- rawdatamatF }
				PLCMPFSLRrawPB[[paste0(sample,'_L_clust_', cluster)]] <- rawPB
			}
		}
	}
	for(sample in unique(PLCRmeta$SampleType)){
		subsetted <- subset(PLCRMPFS, subset=SampleType==sample)
		rawdata <- subsetted$Spatial@data
		rawdatamat <- as.matrix(rawdata)
		if(cluster==3){ 
			cells <- NULL
		} else if(cluster==4){
			cells <- subsetted@meta.data[[param]]==cluster | subsetted@meta.data[[param]]==cluster-1
		} else {
			cells <- subsetted@meta.data[[param]]==cluster
		}
		if(!is.null(cells)){
			if(sum(cells)>1){
				rawdatamatF <- rawdatamat[, cells]
				if(sum(cells)>1){ rawPB <- rowSums(rawdatamatF)} 
				else { rawPB <- rawdatamatF }
				PLCMPFSLRrawPB[[paste0(sample,'_R_clust_', cluster)]] <- rawPB
			}
		}
	}
}

#############
#############

finalcluststouse <- as.numeric(unique(substring(names(ILCMPFSLRrawPB), nchar(names(ILCMPFSLRrawPB)))))

#### ILC
dataframe <- as.data.frame(ILCMPFSLRrawPB)
for(cluster in finalcluststouse){
	print(cluster)
	touse <- grepl(paste0('clust_', cluster), names(dataframe))
	
	df <- dataframe[, touse]
	df2 <- as.data.frame(apply(df, 1:2, round))
	df3 <- df2[rowSums(df2)>=ncol(df2), ]
	
	names <- names(df3)
	underscore <- as.data.frame(gregexpr('_', names))
	first <- c(); second <- c();
	for(i in 1:length(underscore)){
		first <- c(first, underscore[[i]][1])
		second <- c(second, underscore[[i]][2]);
	}
	trmt <- substring(names, first+1, second-1)
	rep <- substring(names, 4, 4)
	side <- substring(names, second+1, second+1)
	lib <- data.frame('Trmt'=trmt, 'Rep'=rep, 'Side'=side, row.names=names)
				
	########## Pro vs Mal
	dds <- DESeqDataSetFromMatrix(countData = df3, colData = lib, design = ~Trmt) #+Rep)
	print(dds)
	dds$Trmt <- relevel(dds$Trmt, ref = 'Mal')  ## For Pro vs Mal, set ref as Mal
	ddsDESeq <- DESeq(dds, test='Wald')
		
	masterres <- as.data.frame(mcols(ddsDESeq))
	masterres <- masterres[,c(1,2,4,7)]	    ## saves the baseMean, baseVar, dispGeneEst and dispersion     
	
	resraw <- results(ddsDESeq, name=resultsNames(ddsDESeq)[3]) ## 3 is the Pro vs Mal comparison
		
	res <- as.data.frame(resraw)
	res <- res[, c(2,3,5,6)]
	res[is.na(res[,1]), 1] <- 0
	res[is.na(res[,3]), 3] <- 1
	res[is.na(res[,4]), 4] <- 1
	names(res) <- gsub('log2FoldChange', 'log2FC', names(res))
	names(res) <- paste0('Pro', '_', names(res))
	
	masterres <- cbind(masterres, res)
	masterresF <- masterres[order(masterres$Pro_pvalue), ]
	
	write.table(masterresF, paste0('DESeq2_ILC_Pro_vs_Mal_Clust_', cluster, '_res.txt'), quote=F, sep='\t')
	
	########## FT vs Sham
	dds <- DESeqDataSetFromMatrix(countData = df3, colData = lib, design = ~Trmt) #+Rep)
	print(dds)
	dds$Trmt <- relevel(dds$Trmt, ref = 'Sham')
	ddsDESeq <- DESeq(dds, test='Wald')
		
	masterres <- as.data.frame(mcols(ddsDESeq))
	masterres <- masterres[,c(1,2,4,7)]	
	
	resraw <- results(ddsDESeq, name=resultsNames(ddsDESeq)[2]) ## 2 is the FT vs Sham comparison
		
	res <- as.data.frame(resraw)
	res <- res[, c(2,3,5,6)]
	res[is.na(res[,1]), 1] <- 0
	res[is.na(res[,3]), 3] <- 1
	res[is.na(res[,4]), 4] <- 1
	names(res) <- gsub('log2FoldChange', 'log2FC', names(res))
	names(res) <- paste0('Pro', '_', names(res))
	
	masterres <- cbind(masterres, res)
	masterresF <- masterres[order(masterres$Pro_pvalue), ]
	
	write.table(masterresF, paste0('DESeq2_ILC_FT_vs_Sham_Clust_', cluster, '_res.txt'), quote=F, sep='\t')
}	

################

#### PLC
dataframe <- as.data.frame(PLCMPFSLRrawPB)
for(cluster in finalcluststouse){
	print(cluster)
	touse <- grepl(paste0('clust_', cluster), names(dataframe))
	
	df <- dataframe[, touse]
	df2 <- as.data.frame(apply(df, 1:2, round))
	df3 <- df2[rowSums(df2)>=ncol(df2), ]
	
	names <- names(df3)
	underscore <- as.data.frame(gregexpr('_', names))
	first <- c(); second <- c();
	for(i in 1:length(underscore)){
		first <- c(first, underscore[[i]][1])
		second <- c(second, underscore[[i]][2]);
	}
	trmt <- substring(names, first+1, second-1)
	rep <- substring(names, 4, 4)
	side <- substring(names, second+1, second+1)
	lib <- data.frame('Trmt'=trmt, 'Rep'=rep, 'Side'=side, row.names=names)
				
	########## Pro vs Mal
	dds <- DESeqDataSetFromMatrix(countData = df3, colData = lib, design = ~Trmt) #+Rep)
	print(dds)
	dds$Trmt <- relevel(dds$Trmt, ref = 'Mal')  ## For Pro vs Mal, set ref as Mal
	ddsDESeq <- DESeq(dds, test='Wald')
		
	masterres <- as.data.frame(mcols(ddsDESeq))
	masterres <- masterres[,c(1,2,4,7)]	    ## saves the baseMean, baseVar, dispGeneEst and dispersion     
	
	resraw <- results(ddsDESeq, name=resultsNames(ddsDESeq)[3]) ## 3 is the Pro vs Mal comparison
		
	res <- as.data.frame(resraw)
	res <- res[, c(2,3,5,6)]
	res[is.na(res[,1]), 1] <- 0
	res[is.na(res[,3]), 3] <- 1
	res[is.na(res[,4]), 4] <- 1
	names(res) <- gsub('log2FoldChange', 'log2FC', names(res))
	names(res) <- paste0('Pro', '_', names(res))
	
	masterres <- cbind(masterres, res)
	masterresF <- masterres[order(masterres$Pro_pvalue), ]
	
	write.table(masterresF, paste0('DESeq2_PLC_Pro_vs_Mal_Clust_', cluster, '_res.txt'), quote=F, sep='\t')
	
	########## FT vs Sham
	dds <- DESeqDataSetFromMatrix(countData = df3, colData = lib, design = ~Trmt) #+Rep)
	print(dds)
	dds$Trmt <- relevel(dds$Trmt, ref = 'Sham')
	ddsDESeq <- DESeq(dds, test='Wald')
		
	masterres <- as.data.frame(mcols(ddsDESeq))
	masterres <- masterres[,c(1,2,4,7)]	
	
	resraw <- results(ddsDESeq, name=resultsNames(ddsDESeq)[2]) ## 2 is the FT vs Sham comparison
		
	res <- as.data.frame(resraw)
	res <- res[, c(2,3,5,6)]
	res[is.na(res[,1]), 1] <- 0
	res[is.na(res[,3]), 3] <- 1
	res[is.na(res[,4]), 4] <- 1
	names(res) <- gsub('log2FoldChange', 'log2FC', names(res))
	names(res) <- paste0('Pro', '_', names(res))
	
	masterres <- cbind(masterres, res)
	masterresF <- masterres[order(masterres$Pro_pvalue), ]
	
	write.table(masterresF, paste0('DESeq2_PLC_FT_vs_Sham_Clust_', cluster, '_res.txt'), quote=F, sep='\t')
}	