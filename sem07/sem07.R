library(VennDiagram)

# Prep

dat <- read.table("bottomly_count_table.tsv", header = TRUE, 
                  row.names = 1)
des <- read.table("bottomly_phenodata.tsv", header = TRUE, 
                  row.names = 1)
all(rownames(des) == colnames(dat))
table(des$strain)

group <- factor(c(rep("1", 10), rep("2", 11)))

# GLM edgeR

library(edgeR)
dge.glm <- DGEList(counts = dat, group = group)
str(dge.glm)
design <- model.matrix(~group)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))

lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

# Mini exercise

# Prep 1
keeps <- apply(dat,1,function(x) {sum(x) > 0})
datFilt <- dat[keeps,]

# Analysis 1
dge.glm <- DGEList(counts = datFilt, group = group)
str(dge.glm)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))

lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

# Prep 2
keeps2 <- apply(dat,1,function(x) {sum(x) > 0})
datFilt2 <- dat[keeps,]

# Analysis 2
dge.glm <- DGEList(counts = datFilt2, group = group)
str(dge.glm)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)

plotBCV(dge.glm.tag.disp)

fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))

lrt <- glmLRT(fit, coef = 2)
topTags(lrt)

tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])

interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))

tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

# DESeq

library(DESeq)
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))

deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)

deSeqDat <- estimateDispersions(deSeqDat)

# plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)

results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)

plotMA(results)

# Voom & limma

library(limma)
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)

dat.voomed

fit <- lmFit(dat.voomed, design)
fit <- eBayes(fit)
topTable(fit)

# Take Home Problem

# Redo voom + limma analysis, specifying the appropriate coef
voomlim <- topTable(fit,coef=2,number=Inf)

# Get results at threshold
thresh <- 1e-10
edgeR_res <- rownames(tt.glm$table[tt.glm$table$FDR < thresh, ])
DEseq_res <- na.omit(results[results$padj < thresh, ])$id
voomlim_res <- rownames(voomlim[voomlim$adj.P.Val < thresh, ])

all_genes <- list(EdgeR=edgeR_res,DESeq=DEseq_res,VoomLim=voomlim_res)
plot.new()
venn_plot <- venn.diagram(all_genes, filename = NULL, fill = c("red", "blue","green"), force.unique=TRUE,ext.text=FALSE,cat.cex=1.5)
grid.draw(venn_plot)