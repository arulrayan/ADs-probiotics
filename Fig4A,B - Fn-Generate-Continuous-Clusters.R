generatecontinuousclusters2 <- function(together = dataframe) {
	#print(str(together))
	### THIS FUNCTION MUST HAVE THE FIRST COLUMN OF THE DATAFRAME ENTITLED 'MODULECOLORS'.
	distances <- function(x, y){
		vec = NULL
		for(i in 1:ncol(y)){
			vecy <- y[,i]
			vec <- rbind(vec, stats::cor(x, vecy))
		}
		return(vec)
	}
	
	meanofmodules <- NULL
	modules <- unique(together$modulecolors)
	for(module in modules){
		temp <- together[together$modulecolors==module, ]
		temp <- temp[, 2:ncol(temp)]	
		average <- unlist(lapply(temp, mean))
		meanofmodules <- rbind(meanofmodules, average)
	}
	rownames(meanofmodules) <- modules
	
	meanofmodules <- scale(meanofmodules)
	
	orderofmodules <- NULL
	namesofmodules <- NULL

	moduleone <- 1
	moduletempt <- as.data.frame(t(meanofmodules))
	
	for(i in 1:(ncol(moduletempt)-1)){
		orderofmodules <- c(orderofmodules, moduleone)
		namesofmodules <- c(namesofmodules, names(moduletempt)[moduleone])
		seq <- 1:ncol(moduletempt)
		if(length(seq)>2){
			newseq <- seq[seq!=moduleone]
			moduleonevec <- moduletempt[, moduleone]
			moduletempt <- moduletempt[, newseq]
			
			disstats <- distances(moduleonevec, moduletempt)
			rownames(disstats) <- names(moduletempt) 
			moduleone <- which(disstats==max(disstats))
			#corrstats <- stats::cor(moduletempt, moduleonevec, method='spearman')
		} else{
			newseq <- seq[seq!=moduleone]
			namesofmodules <- c(namesofmodules, names(moduletempt)[newseq])
			print(namesofmodules)
		}
	}

	together2 <- NULL
	for(module in namesofmodules){
		temp <- together[together$modulecolors==module, ]
		together2<- rbind(together2, temp)	
	}
	
	Annotationrow2 <- as.character(together2$modulecolors)
	Annotationrow2 <- as.data.frame(Annotationrow2, row.names = rownames(together2), stringsAsFactors=F)
	names(Annotationrow2) <- 'module2'
	
	newlog2FC <- together2[, 2:(ncol(together2))]
	names(namesofmodules) <- namesofmodules
	
	listtoreturn <- list()
	listtoreturn[['modulemeans']] <- t(meanofmodules)
	listtoreturn[['together2']] <- together2
	listtoreturn[['Annotationrow2']] <- Annotationrow2
	listtoreturn[['modulenames']] <- namesofmodules
	
	return(listtoreturn)
}
