library(anRichment)

###################

setwd('/path/to/my/directory/')

allresults <- as.data.frame(read.table(file='All DESeq2 w Rep APEGLM-ILCsLC, with Pro-Mal, FC=40%, sym, padj=0,01.txt', header=TRUE, sep='\t'))

log2FC <- lapply(names(allresults), function(ch) grep('log2FC', ch))
log2FC <- allresults[, which(log2FC==1)] 

allabssums <- allresults[, 290]

###################

log2FCsubset <- log2FC[allabssums>=6, ]
nclusters <- 18

log2FCscaled <- as.data.frame(scale(log2FCsubset))
clusternames <- c(1:nclusters)
for(m in 1:nclusters) { clusternames[m] <- paste('Cluster-', m, sep="") }
	
set.seed(1234)
km.res <- kmeans(log2FCscaled, nclusters, iter.max=300, nstart = 300, algorithm = "Lloyd")

centers <- t(km.res$centers)
clusters <- km.res$cluster
clusters2 <- data.frame('Cluster number' = paste0('Cluster ', clusters), 'genenames' =attr(clusters, 'names'))

write.table(clusters2, file='Cluster identities.txt', quote=F, sep='\t')
clustersordered <- clusters2[order(clusters2$Cluster.number), ]
write.table(clustersordered, file='Cluster identities ordered by clust.txt', quote=F, row.names=F, sep='\t')

###################  ANRICHMENT

genenames <- clustersordered$genenames
finalclusters <- clustersordered$Cluster.number 
GOcollection = buildGOcollection(organism = "rat")
biosysCollection = BioSystemsCollection("rat")
internalColl = internalCollection(organism = "rat");
combinedCollection = mergeCollections(GOcollection, internalColl, biosysCollection)

genenamesentrez <- convert2entrez(organism = 'rat', symbol = genenames)

GOenrichment = enrichmentAnalysis(classLabels = finalclusters, identifiers = genenamesentrez, refCollection = combinedCollection, useBackground = "given", threshold = 1e-4, thresholdType = "Bonferroni", getOverlapEntrez = F, getOverlapSymbols = TRUE, maxReportedOverlapGenes=300, geneSeparator='/', ignoreLabels = "grey");

enrichmentresults <- GOenrichment$enrichmentTable
enrichmentresults <- enrichmentresults[order(enrichmentresults$FDR), ]
write.csv(enrichmentresults, file = paste(nclusters, "_clust-KM-mod-anRich-combinedcoll.csv", sep=''), row.names = FALSE)		
