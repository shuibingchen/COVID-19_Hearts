# use edgeR package to analysis the gene reads and generate MDS plot and sample distance heatmap


# load the counts of each sample
MOCK<-read.table(file = "mock_gene_count.tsv",header = TRUE,row.names = 1)

IN<-read.table(file = "infect_gene_count.tsv",header = TRUE,row.names = 1)


MOCK_IN<-cbind(MOCK,IN)

# rename the colname to the group name
group<-rep(c("mock","SARSCOV2"),each=3)

group<-as.factor(group)

colnames(MOCK_IN)<-paste(group,1:3,sep = "-")

# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(MOCK_IN)

y<-DGEList(counts=MOCK_IN,genes=genelist)

# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]


# Normalization
y <- calcNormFactors(y, method="TMM")

y$samples$group <- group

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# exploring differences between libraries
# MDS plot
pch<-c(15,16)
colors<-rep(c("red","green"),3)

plotMDS(y,col=colors[group],pch=pch[group])

legend("top",legend=levels(group),pch=pch,col=colors)


# samples distance plot
library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(MOCK_IN))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(MOCK_IN))
colnames(sampleDistMatrix) <- paste(colnames(MOCK_IN))
colors2 <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sample_distance_plot<-
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors2)

