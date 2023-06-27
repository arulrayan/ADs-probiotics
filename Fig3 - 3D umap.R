library(plotly)
library(umap)
library(processx)

setwd('/path/to/my/directory/')

allresults <- as.data.frame(read.table(file='All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt', header=TRUE, sep='\t'))

names(allresults) <- gsub('LConly', 'LC', names(allresults))

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

allabssums <- lapply(names(allresults), function(ch) grep('allabssums', ch))
allabssums <- allresults[, which(allabssums==1)]

################

colors <- list(regions = c(BLA = "#999999", CGC="#000000", dorDG="#E69F00", ILC="#56B4E9", LC="#009E73", mPOA="green", NACShell="#F0E442", PLC="#0072B2", Raphe="#D55E00", venDG="#CC79A7"))

regioncolors <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "green", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") 
trmtchoices <- c('circle', 'square', 'diamond', 'x', 'cross')
sizechoices <- c(Bup=21, Des=21, FT=18, Mal=7, Pro=18)

## UMAP parameters
neighbors=6 
mindist=0.1 

##PCA components
numcomps=12
centered = 1
scaled = 1

## Number of region-treatment pairs a gene must be differentially expressed in to be included in plot
X=6

################

log2FCF <- log2FC[allabssums>=X, ]	
pca <- prcomp(t(log2FCF), center=as.logical(centered), scale=as.logical(scaled))
PCAs <- pca$x[, 1:numcomps]

PCAsUmap <- umap(PCAs, random_state=123, n_components=3, n_neighbors=neighbors, min_dist=mindist)
PCAslayout <- as.data.frame(PCAsUmap$layout)

################
### PLOTTING MAIN FIGURE (FIG. 2A)

rownames(PCAslayout) <- gsub('_log2FC', '', rownames(PCAslayout))

samples <- rownames(PCAslayout)
underscore <- unlist(gregexpr('_', samples))
regions <- substring(samples, 1, underscore-1)
treatments <- substring(samples, underscore+1, nchar(samples))
PCAslayout <- cbind(PCAslayout, regions, treatments)

sizestouse <- unlist(sizechoices[match(treatments, names(sizechoices))])
PCAslayout <- cbind(PCAslayout, sizestouse)
names(PCAslayout)[1:3] <- c('UMAP1', 'UMAP2', 'UMAP3')

f1 <- list(family = "Arial, sans-serif", size = 20, color = "black")
f2 <- list(family = "Arial, sans-serif", size = 14, color = "black")
ax <- list(title = "UMAP1", titlefont = f1, tickfont = f2, gridcolor = toRGB("gray50"), gridwidth = 1, linecolor = toRGB("black"), linewidth = 4)
ay <- list(title = "UMAP2", titlefont = f1, tickfont = f2, gridcolor = toRGB("gray50"), gridwidth = 1, linecolor = toRGB("black"), linewidth = 4)
az <- list(title = "UMAP3", titlefont = f1, tickfont = f2, gridcolor = toRGB("gray50"), gridwidth = 1, linecolor = toRGB("black"), linewidth = 4 )

fig <- plot_ly(PCAslayout, x=~UMAP1, y=~UMAP2, z=~UMAP3, type='scatter3d', mode='markers', stroke=10, strokes='black', color=~regions, colors=regioncolors, symbol=~treatments, symbols=trmtchoices, size=~sizestouse, sizes=c(200,1500))
fig <- fig %>% layout(scene = list(xaxis=ax, yaxis=ay, zaxis=az))

################
### PLOTTING FIG. 2B-D
### CGC mPOA PLC Raphe venDG
regioncolors1 <- c("#F0F0F0", "#000000", "#F0F0F0", "#F0F0F0", "#F0F0F0", "green", "#F0F0F0", "#0072B2", "#D55E00", "#CC79A7")

fig1 <- plot_ly(PCAslayout, x=~UMAP1, y=~UMAP2, z=~UMAP3, type='scatter3d', mode='markers', stroke=10, strokes='black', color=~regions, colors=regioncolors1, symbol=~treatments, symbols=trmtchoices, size=~sizestouse, sizes=c(200,1500)) 
fig1 <- fig1 %>% layout(scene = list(xaxis=ax, yaxis=ay, zaxis=az))

################
### BLA, ILC, LC
regioncolors2 <- c("#999999", "#F0F0F0", "#F0F0F0", "#56B4E9", "#009E73", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0F0F0") 

fig2 <- plot_ly(PCAslayout, x=~UMAP1, y=~UMAP2, z=~UMAP3, type='scatter3d', mode='markers', stroke=10, strokes='black', color=~regions, colors=regioncolors2, symbol=~treatments, symbols=trmtchoices, size=~sizestouse, sizes=c(200,1500)) 
fig2 <- fig2 %>% layout(scene = list(xaxis=ax, yaxis=ay, zaxis=az))

################
### NACshell, dorDG
regioncolors3 <- c("#F0F0F0", "#F0F0F0", "#E69F00", "#F0F0F0", "#F0F0F0", "#F0F0F0", "#F0E442", "#F0F0F0", "#F0F0F0", "#F0F0F0") 

fig3 <- plot_ly(PCAslayout, x=~UMAP1, y=~UMAP2, z=~UMAP3, type='scatter3d', mode='markers', stroke=10, strokes='black', color=~regions, colors=regioncolors3, symbol=~treatments, symbols=trmtchoices, size=~sizestouse, sizes=c(200,1500)) 
fig3 <- fig3 %>% layout(scene = list(xaxis=ax, yaxis=ay, zaxis=az))
