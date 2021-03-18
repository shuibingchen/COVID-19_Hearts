# run_DESeq2.R
# perform differential expression analysis on the 8 bulkRNAseq samples using DESeq2
# 

library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(BiocParallel)
library(apeglm)

workdir <- "./"
countsfile <- file.path(workdir, "source", "raw_counts_genes.clean.txt")
figdir <- file.path(workdir, "figures")
numdir <- file.path(workdir, "info")

# false discovery rate cutoff
alpha <- 0.1

# read in raw counts file
countData <- read.table(file=countsfile, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep='\t', row.names=1)

# create experimental label
colData <- data.frame(condition=factor(c(rep("Healthy",5),rep("Patient",3))))
rownames(colData) <- colnames(countData)

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
			# moderated log2 fold changes
			resLFC <- lfcShrink(dds, res=res, type="ashr", parallel=TRUE, BPPARAM=MulticoreParam(4))
			# reorder results by adjusted-pvalue
			resLFCOrdered <- resLFC[order(resLFC$padj),]
			# make MA-plot
			png(file=file.path(figdir, paste("MA-plot", allconditions[i], "vs", allconditions[j], "LFC", "png", sep=".")), width = 1920, height = 1920, units = "px", pointsize = 6, res=600)
			plotMA(resLFC, main=paste(allconditions[i], "vs", allconditions[j], sep=" "))
			abline(h=c(-1,1), col="dodgerblue", lwd=2)
			dev.off()
			# exporting results to file
			write.table(as.data.frame(resLFCOrdered), file=paste(numdir, file.path("DE", allconditions[i], "vs", allconditions[j], "LFC", "tsv", sep=".")), quote=FALSE, sep='\t', col.names=NA)
		}
	}
}

# data visualization
rld <- rlog(dds, blind=F)

# heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
png(file=file.path(figdir, "heatmap.sample-to-sample.rld.png"), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=hmcol, fontsize_row=6, fontsize_col=6)
dev.off()

# MDS (multidimensional scaling) plot
mds <- data.frame(cmdscale(sampleDistMatrix))
colData <- data.frame(condition=factor(c(rep("Healthy",5),rep("Patient",3))))
rownames(colData) <- rownames(mds)
mds <- cbind(mds, colData)
mds$name <- rownames(mds)

g <- ggplot(mds, aes(x=X1, y=X2, color=condition, label=name))
g <- g + geom_point(shape=19, size=2)
g <- g + xlab("coordinate 1") + ylab("coordinate 2")
g <- g + theme_classic()
ggsave(file=file.path(figdir, "MDS.rld.png"), width=6, height=4.5, dpi=600)

# save R object
saveRDS(dds, file=file.path(numdir, "dds.rds"))

sessionInfo()

