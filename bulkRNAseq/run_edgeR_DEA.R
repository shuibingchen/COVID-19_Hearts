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

# the DEA result for all the genes
# dea <- lrt$table

y <- estimateDisp(y, design, robust = TRUE)

fit<-glmQLFit(y,design,robust = TRUE)

IN_vs_MOCK<-makeContrasts(SARSCOV2-mock,levels = design)

res<-glmQLFTest(fit,contrast = IN_vs_MOCK)

toptag <- topTags(res, n = nrow(y$genes), p.value = 1)

dea <- toptag$table 

dea <- dea[order(dea$FDR, -abs(dea$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

write.table(dea,file = "MOCK_IN_DEA.tsv",quote=FALSE,sep = "\t")
