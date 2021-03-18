# run_DESeq2.R
# perform DESeq differential analysis on the 16 samples (4 group)
# 

library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(BiocParallel)
library(apeglm)

workdir <- "/data/gc-core/taz2008/RNAseq/Project_Chen-XT-8886_200528_ANALYSIS/DESeq2.patient"
countsfile <- paste(workdir, "source", "raw_counts_genes.clean.txt", sep="/")
figdir <- paste(workdir, "figures", sep="/")
numdir <- paste(workdir, "info", sep="/")

# false discovery rate cutoff
alpha <- 0.1

# read in raw counts file
countData <- read.table(file=countsfile, header=TRUE, check.names=FALSE)
rownames(countData) <- countData[,1]
countData <- countData[,-1]
head(countData,n=2)

# create experimental label
colData <- data.frame(condition=factor(c(rep("Healthy",5),rep("Patient",3))))
rownames(colData) <- colnames(countData)
colData

# use parallelization
register(MulticoreParam(4))

# construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
dds

# update factor levels
dds$condition <- factor(dds$condition, levels=c("Healthy","Patient"))

# run DESeq
dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(4))

allconditions <- c("Patient","Healthy")
for (i in c(1:length(allconditions))){
	for (j in c(1:length(allconditions))){
		if (j>i){
			print(paste(allconditions[i], allconditions[j], sep=" v.s. "))
			# DE
			res <- results(dds, contrast=c("condition", allconditions[i], allconditions[j]), alpha=alpha, parallel=TRUE, BPPARAM=MulticoreParam(4))
			# reorder results by adjusted-pvalue
			resOrdered <- res[order(res$padj),]
			summary(res)
			# how many ajusted p-values were less than 0.05?
			sum(res$padj<0.05, na.rm=TRUE)
			# get information about which variables and tests were used
			mcols(res, use.names=TRUE)$description
			# make MA-plot
			png(file=paste(figdir, paste("MA-plot", allconditions[i], "vs", allconditions[j], "png", sep="."), sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
			plotMA(res, main=paste(allconditions[i], "vs", allconditions[j], sep=" "))
			abline(h=c(-1,1), col="dodgerblue", lwd=2)
			dev.off()
			# exporting results to file
			write.table(as.data.frame(resOrdered), file=paste(numdir, paste("DE", allconditions[i], "vs", allconditions[j], "tsv", sep="."), sep="/"), quote=FALSE, sep='\t', col.names=NA)
			# select significantly changed genes only
			##resSig <- subset(resOrdered, padj<0.05)
			##as.data.frame(resSig)
			# independent filtering
			##attr(res,"filterThreshold")
			#metadata(res)$filterNumRej
			png(file=paste(figdir, paste("independent_filtering", allconditions[i], "vs", allconditions[j], "png", sep="."), sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
			##plot(attr(res, "filterNumRej"), type= "b", ylab="number of rejections", xlab="quantiles of filter")
			#lines(attr(res, "lo.fit"), col="red")
			#abline(v=attr(res, "filterTheta"))
			plot(metadata(res)$filterNumRej, type="b", ylab="number of rejections", xlab="quantiles of filter")
			lines(metadata(res)$lo.fit, col="red")
			abline(v=metadata(res)$filterTheta)
			dev.off()
			# test log2FoldChange above or below 1.0
			##resGA <- results(dds, contrast=c("condition", allconditions[i], allconditions[j]), lfcThreshold=1, altHypothesis="greaterAbs")
			##summary(resGA)
			##resGAOrdered <- resGA[order(resGA$padj),]
			##write.table(as.data.frame(resGAOrdered), file=paste(numdir, paste("DE", "fc1.0", allconditions[i], "vs", allconditions[j], "tsv", sep="."), sep="/"), quote=FALSE, sep='\t', col.names=NA)
			# moderated log2 fold changes
			#print(resultsNames(dds))
			#resLFC <- lfcShrink(dds, coef=paste("condition",allconditions[i],"vs",allconditions[j],sep="_"), type="apeglm", parallel=TRUE, BPPARAM=MulticoreParam(4))
			resLFC <- lfcShrink(dds, res=res, type="ashr", parallel=TRUE, BPPARAM=MulticoreParam(4))
			# reorder results by adjusted-pvalue
			resLFCOrdered <- resLFC[order(resLFC$padj),]
			#summary(res)# should be same as above, since p-values are not changed
			# how many ajusted p-values were less than 0.05?
			#sum(resLFC$padj<0.05, na.rm=TRUE)# should be same as above
			# make MA-plot
			png(file=paste(figdir, paste("MA-plot", allconditions[i], "vs", allconditions[j], "LFC", "png", sep="."), sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
			plotMA(resLFC, main=paste(allconditions[i], "vs", allconditions[j], sep=" "))
			abline(h=c(-1,1), col="dodgerblue", lwd=2)
			dev.off()
			# exporting results to file
			write.table(as.data.frame(resLFCOrdered), file=paste(numdir, paste("DE", allconditions[i], "vs", allconditions[j], "LFC", "tsv", sep="."), sep="/"), quote=FALSE, sep='\t', col.names=NA)
		}
	}
}

# plot counts for a given gene
#d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
###ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1, h=0)) + scale_y_log10(breaks=c(25,100,400))
#ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1, h=0)) + scale_y_log10()

# data visualization
rld <- rlog(dds, blind=F)
vsd <- vst(dds, blind=F)
head(assay(rld), n=3)
head(assay(vsd), n=3)

# effect of transformation
notAllZero <- (rowSums(counts(dds))>0)
png(file=paste(figdir, "meanSdPlot.log.png", sep="/"), width = 2400, height = 1920, units = "px", pointsize = 6, res=600)
meanSdPlot(log2(counts(dds, normalized=TRUE)[notAllZero,]+1))
dev.off()
png(file=paste(figdir, "meanSdPlot.rld.png", sep="/"), width = 2400, height = 1920, units = "px", pointsize = 6, res=600)
meanSdPlot(assay(rld[notAllZero,]))
dev.off()
png(file=paste(figdir, "meanSdPlot.vsd.png", sep="/"), width = 2400, height = 1920, units = "px", pointsize = 6, res=600)
meanSdPlot(assay(vsd[notAllZero,]))
dev.off()

# heatmap of count matrix
# for top20 expressed genes over all samples
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:20]
##hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol <- colorRampPalette(brewer.pal(9, "Oranges"))(100)
png(file=paste(figdir, "heatmap.top20.counts.log.png", sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
heatmap.2(log2(counts(dds, normalized=TRUE)[select, ]+1), col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(12,12))
dev.off()
#df <- as.data.frame(colData(dds)[,c("condition")])
png(file=paste(figdir, "heatmap.top20.counts.rld.png", sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
heatmap.2(assay(rld)[select,], col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(12,12))
dev.off()
png(file=paste(figdir, "heatmap.top20.counts.vsd.png", sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
heatmap.2(assay(vsd)[select,], col=hmcol, Rowv=FALSE, Colv=FALSE, scale="none", dendrogram="none", trace="none", margin=c(12,12))
dev.off()

# heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(rld$condition, rld$pair, sep="-")
#colnames(sampleDistMatrix) <- NULL
##hc <- hclust(sampleDists)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
png(file=paste(figdir, "heatmap.sample-to-sample.rld.png", sep="/"), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
##heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace='none', col = rev(hmcol), margin=c(12, 12))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=hmcol, fontsize_row=6, fontsize_col=6)
dev.off()

# PCA
###png(file=paste(figdir, "PCA.rld.png", sep="/"), width = 1920, height = 1280, units = "px", pointsize = 6, res=600)
png(file=paste(figdir, "PCA.rld.png", sep="/"), width = 3600, height = 3200, units = "px", res=600)
plotPCA(rld, intgroup=c("condition"))
dev.off()

# updated PCA plot with sample labels
mypcadata <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
head(mypcadata,n=3)
percentVar <- round(100*attr(mypcadata, "percentVar"))
g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_point(size=3, shape=19) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)
g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_text() + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.labeled.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)
write.table(mypcadata, paste(numdir, "mypcadata.txt", sep="/"), quote=F, row.names = T, col.names = NA, sep="\t")

# only use virus genes for PCA
virus.genes <- grep("^GU280", rownames(rld), value=T)
mypcadata <- plotPCA(rld[virus.genes,], intgroup=c("condition"), returnData=TRUE)
head(mypcadata,n=3)
percentVar <- round(100*attr(mypcadata, "percentVar"))
g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_point(size=3, shape=19) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.virus_genes_only.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)
g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
g <- g + geom_text() + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
g <- g + theme_classic()
ggsave(file=paste(figdir, "PCA.rld.v2.virus_genes_only.labeled.png", sep="/"), width = 8, height = 6, type = "cairo", dpi = 600)

# dispersion plot
png(file=paste(figdir, "dispersion.png", sep="/"), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
plotDispEsts(dds)
dev.off()

# boxplot of Cook's distances
png(file=paste(figdir, "cooks.boxplot.png", sep="/"), width = 2400, height = 1920, units = "px", pointsize = 6, res=600)
par(mar=c(10,4,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
title("Cook's distances")
dev.off()

# supplment: MDS (multidimensional scaling) plot
mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, colData(rld))
###png(file=paste(figdir, "MDSplot.rld.png", sep="/"), width = 1920, height = 1280, units = "px", pointsize = 6, res=600)
png(file=paste(figdir, "MDSplot.rld.png", sep="/"), width = 4800, height = 3200, units = "px", res=600)
####qplot(X1, X2, color=condition, shape=pair, data=mds)+scale_shape_manual(values=LETTERS[1:12])+labs(shape="patient")
qplot(X1, X2, color=condition, data=as.data.frame(mds))
dev.off()

# raw counts and normalized counts
rawcounts <- counts(dds)
normcounts <- counts(dds, normalized=TRUE)
wh.rows.norm <- match(rownames(rawcounts), rownames(normcounts))
combResults <- cbind(rawcounts, normcounts[wh.rows.norm, ])
write.table(as.data.frame(combResults), file=paste(numdir, "details_counts.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)

# regularized log transformation
write.table(as.data.frame(assay(rld)), file=paste(numdir, "details_rld.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)
# variance stabilizing transformation
write.table(as.data.frame(assay(vsd)), file=paste(numdir, "details_vsd.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)

# size factors
write.table(as.data.frame(sizeFactors(dds)), file=paste(numdir, "details_size_factors.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)

#!!! NEW
# remove batch effect
#condition <- colData$condition
#design <- model.matrix(~condition)
#rlogMat <- limma::removeBatchEffect(assay(rld), batch=colData$batch, design=design)
#assay(rld) <- rlogMat
#assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch=vsd$batch, design=design)

# updated PCA plot with sample labels
#mypcadata <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
#head(mypcadata,n=3)
#percentVar <- round(100*attr(mypcadata, "percentVar"))
#g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
#g <- g + geom_point(size=3, shape=19) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
#g <- g + theme_classic()
#ggsave(file=paste(figdir, "PCA.rld.no_batch_effect.v2.png", sep="/"), width = 6, height = 6, type = "cairo", dpi = 600)
#g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
#g <- g + geom_text() + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
#g <- g + coord_fixed()
#g <- g + theme_classic()
#ggsave(file=paste(figdir, "PCA.rld.no_batch_effect.v2.labeled.png", sep="/"), width = 6, height = 6, type = "cairo", dpi = 600)

# heatmap of the sample-to-sample distances
#sampleDists <- dist(t(assay(rld)))
#sampleDistMatrix <- as.matrix(sampleDists)
##rownames(sampleDistMatrix) <- paste(rld$condition, rld$pair, sep="-")
##colnames(sampleDistMatrix) <- NULL
##hc <- hclust(sampleDists)
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
#png(file=paste(figdir, "heatmap.sample-to-sample.rld.no_batch_effect.png", sep="/"), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
##heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc), symm=TRUE, trace='none', col = rev(hmcol), margin=c(12, 12))
#pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=hmcol, fontsize_row=6, fontsize_col=6)
#dev.off()

# supplment: MDS (multidimensional scaling) plot
#mds <- data.frame(cmdscale(sampleDistMatrix))
#mds <- cbind(mds, colData(rld))
###png(file=paste(figdir, "MDSplot.rld.png", sep="/"), width = 1920, height = 1280, units = "px", pointsize = 6, res=600)
#png(file=paste(figdir, "MDSplot.rld.no_batch_effect.png", sep="/"), width = 4800, height = 3200, units = "px", res=600)
####qplot(X1, X2, color=condition, shape=pair, data=mds)+scale_shape_manual(values=LETTERS[1:12])+labs(shape="patient")
#qplot(X1, X2, color=condition, data=as.data.frame(mds))
#dev.off()

# regularized log transformation
#write.table(as.data.frame(assay(rld)), file=paste(numdir, "details_rld.no_batch_effect.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)
# variance stabilizing transformation
#write.table(as.data.frame(assay(vsd)), file=paste(numdir, "details_vsd.no_batch_effect.tsv", sep="/"), quote=FALSE, sep='\t', col.names=NA)

# save R object
saveRDS(dds, file=paste(numdir, "dds.rds", sep="/"))
saveRDS(rld, file=paste(numdir, "rld.rds", sep="/"))
saveRDS(vsd, file=paste(numdir, "vsd.rds", sep="/"))

sessionInfo()

