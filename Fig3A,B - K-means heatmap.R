library(Seurat)
library(pheatmap)
source('/path/to/my/directory/Fig3A,B - Fn-Generate-Continuous-Clusters.R')

generateheatmap <- function(x, colorscheme, Annotationnames, Annotationrow, colors, breaks, clustercols, clusterrows){
	if(nrow(x)>=5000){
		cellheight=0.09
	} else if(nrow(x)>=2500){
		cellheight=0.15
	} else if(nrow(x)>=1250){
		cellheight=0.225
	} else { cellheight=0.3 }
	heatmap <- pheatmap(x, color=colorscheme, annotation_col=Annotationnames, annotation_row=Annotationrow, annotation_colors=colors, breaks=breaks, cluster_cols=as.logical(clustercols), cluster_rows=as.logical(clusterrows), clustering_distance_cols='correlation', show_rownames=FALSE, border_color='NA', cellwidth=15, cellheight=cellheight, fontsize=14, angle_col=90, annotation_names_col=F)
	return(heatmap)
}
save_pheatmap_png <- function(x, filename, y, totalcolors, res = 300) {
	if(nrow(y)>=5000){
		cellheight=0.09
	} else if(nrow(y)>=2500){
		cellheight=0.15
	} else if(nrow(y)>=1250){
		cellheight=0.225
	} else { cellheight=0.3 }
	imageheight=max(cellheight*length(row.names(y))+4000, 150*totalcolors)
    png(filename, width = 80*length(names(y))+1000, height = imageheight, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}
breaksz = seq(-1, 1, 2/100)
Seuratcolor <- PurpleAndYellow(100)  ### generates Seurat Magenta-Black-Yellow scheme

#####################

setwd('/path/to/my/directory/')

allresults <- as.data.frame(read.table(file=file.path(currdir, 'All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt'), header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

baseMean <- lapply(names(allresults), function(ch) grep('baseMean', ch))
baseMean <- allresults[, which(baseMean==1)] 

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

allabssums <- allresults[, 290]

##################### GETTING THE SUBSETTED MATRIX FOR HEATMAP PLOTTING

data <- read.table(file='Cluster identities ordered by clust.txt', header=T, sep='\t')
clusters <- data$Cluster.number
genenames <- data$genenames 

log2FCF <- log2FC
log2FCF <- log2FCF[match(genenames, rownames(log2FCF)), ]

##################### PRE-PROCESSING FOR DRAWING LABELS FOR REGIONS AND TREATMENTS IN PHEATMAP

names(log2FCF) <- gsub('_log2FC', "", names(log2FCF))
RTPnames <- unlist(names(log2FCF))
underscore <- unlist(gregexpr('_', RTPnames))
Regions <- substring(RTPnames, 1, underscore-1)
Treatments <- substring(RTPnames, underscore+1)
Annotationnames <- data.frame('Regions'=Regions, 'Treatments'=Treatments, row.names= names(log2FCF))

Annotationrow <- as.data.frame(as.factor(clusters), row.names = rownames(log2FCF), stringsAsFactors=F)
names(Annotationrow) <- 'Module'

modulecolors <- standardColors(length(unique(clusters)))
names(modulecolors) <- as.factor(unique(clusters))

colors <- list(Regions = c(BLA = "#999999", CGC="#000000", dorDG="#E69F00", ILC="#56B4E9", LC="#009E73", mPOA="green", NACShell="#F0E442", PLC="#0072B2", Raphe="#D55E00", venDG="#CC79A7"), Treatments = c(Bup='#A04700', Des='#1C91D4', FT='#007756', Mal='#FFCB57', Pro='#666666'), 'Module' = modulecolors)

##################### GENERATE 'CONTINUOUS'-LOOKING CLUSTERS

together <- cbind(clusters, log2FCF)	
names(together)[1] <- 'modulecolors'
continuousclusters <- generatecontinuousclusters2(together=together)
continuouslog2FC <- continuousclusters$together2
continuouslog2FC <- continuouslog2FC[, 2:ncol(continuouslog2FC)]
Annotationrow2 <- continuousclusters$Annotationrow2
names(Annotationrow2) <- 'Module'

##################### GENERATE HEATMAP TO SEE COLUMN ASSIGNEMENTS

log2FCcontinuousheatmap <- generateheatmap(continuouslog2FC, colorscheme=Seuratcolor, Annotationnames, Annotationrow2, colors, breaksz, clustercols=1, clusterrows=0)

##################### BREAK UP COLUMNS ACCORDING TO UMAP PARTITIONING, FIGURE 2

contcolumnorder <- log2FCcontinuousheatmap$tree_col$order
contcolumn1 <- c(16, 18, 15, 17, 19, 24, 20, 22, 21, 23, 5, 1, 4, 2, 3)
contcolumn2 <- c(9, 28, 38, 43, 48, 6, 25, 34, 39, 44, 7, 26, 36, 41, 46, 35, 40, 45, 8, 27, 37, 42, 47)
contcolumn3 <- c(33, 29, 31, 10, 12, 13, 32, 30, 11, 14)

filename <- 'KM_18_clusts_'
contfilename1 <- file.path(currdir, paste(filename, "cont-sc-phtmp2-SeuV2-Grp1-V2.png", sep='-'))
contfilename2 <- file.path(currdir, paste(filename, "cont-sc-phtmp2-SeuV2-Grp2-V2.png", sep='-'))
contfilename3 <- file.path(currdir, paste(filename, "cont-sc-phtmp2-SeuV2-Grp3-V2.png", sep='-'))

log2FCFcontinuousheatmapgroup1 <- generateheatmap(continuouslog2FC[, contcolumn1], colorscheme=Seuratcolor, Annotationnames, Annotationrow, colors, breaksz, clustercols=0, clusterrows=0)
log2FCFcontinuousheatmapgroup2 <- generateheatmap(continuouslog2FC[, contcolumn2], colorscheme=Seuratcolor, Annotationnames, Annotationrow, colors, breaksz, clustercols=0, clusterrows=0)
log2FCFcontinuousheatmapgroup3 <- generateheatmap(continuouslog2FC[, contcolumn3], colorscheme=Seuratcolor, Annotationnames, Annotationrow, colors, breaksz, clustercols=0, clusterrows=0)
save_pheatmap_png(x=log2FCFcontinuousheatmapgroup1, filename=contfilename1, y=continuouslog2FC[, contcolumn1], totalcolors=length(unlist(colors)))
save_pheatmap_png(x=log2FCFcontinuousheatmapgroup2, filename=contfilename2, y=continuouslog2FC[, contcolumn2], totalcolors=length(unlist(colors)))
save_pheatmap_png(x=log2FCFcontinuousheatmapgroup3, filename=contfilename3, y=continuouslog2FC[, contcolumn3], totalcolors=length(unlist(colors)))	