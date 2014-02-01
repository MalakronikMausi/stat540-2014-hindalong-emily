# Read in data and design file
library(lattice)
prDes <- read.table("../GSE4051_design.tsv", header = TRUE, as.is = 1)
prDat <- read.table("../GSE4051_data.tsv")

# Select a gene at random and bind with design info
set.seed(987)
theGene <- sample(1:nrow(prDat),1)
pDat <- data.frame(prDes, gExp=unlist(prDat[theGene,]))

# Inspect the data
aggregate(gExp ~ gType, pDat, FUN = mean)
stripplot(gType ~ gExp, pDat)

# t-test w/ different variance
ttRes_diff_var <- t.test(gExp ~ gType, pDat)
ttRes_diff_var$statistic
ttRes_diff_var$p.value

# t-test w/ same variance
ttRes_same_var <- t.test(gExp ~ gType, pDat, var.equal = TRUE)
ttRes_same_var$statistic
ttRes_same_var$p.value

# KS
ksRes <- suppressWarnings(ks.test(unlist(subset(pDat,subset=gType=="wt",select = c("gExp"))),unlist(subset(pDat,subset=gType=="NrlKO",select = c("gExp")))))
ksRes$statistic
ksRes$p.value

# Wilcox
wilcoxRes <- suppressWarnings(wilcox.test(gExp ~ gType, pDat))
wilcoxRes$statistic
wilcoxRes$p.value

# Read in mini data
library(plyr)
kDat <- read.table("../GSE4051_MINI.txt",header=TRUE,row.names=1)

# Aggregate
aggregate(eggBomb ~ gType * devStage, kDat, FUN = range)

# plyr
keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at",
               "1416119_at", "1432141_x_at", "1429226_at" )
miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)
miniDat <- data.frame(gExp = as.vector(t(as.matrix(miniDat))),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = keepGenes))
miniDat <- suppressWarnings(data.frame(prDes, miniDat))

ttRes1 <- d_ply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x))
ttRes2 <- dlply(miniDat, ~ gene, function(x) t.test(gExp ~ gType, x))
ttRes2

# ***** Take home work *****

# p-values for 100 genes for three different test statistics

samp100 <- sample(1:nrow(prDat),100)
subNames <- rownames(prDat)[samp100]
subDat <- prDat[samp100,]
subDat <- data.frame(gExp = as.vector(t(as.matrix(subDat))),
                      gene = factor(rep(rownames(subDat), each = ncol(subDat)),
                                    levels = subNames))
subDat <- suppressWarnings(data.frame(prDes, subDat))

ttRes <- ddply(subDat, ~ gene, function(z) {
  zz <- t.test(gExp ~ gType, z)
  zzz <- suppressWarnings(wilcox.test(gExp ~ gType, z))
  zzzz <- suppressWarnings(ks.test(unlist(subset(z,subset=gType=="wt",select = c("gExp"))),unlist(subset(z,subset=gType=="NrlKO",select = c("gExp")))))
  round(c(pVal_t = zz$p.value, pVal_wilcox = zzz$p.value, pVal_ks = zzzz$p.value), 4)
})
ttRes

# Plot the relationship - ks seems bit off from the others, I guess because it is sensitive to the shape of the distribution?
xyplot(pVal_wilcox + pVal_ks ~ pVal_t,ttRes,outer=TRUE)
xyplot(pVal_ks ~ pVal_wilcox,ttRes)

# Log transform doesn't help much
xyplot(pVal_wilcox + pVal_ks ~ pVal_t,ttRes,outer=TRUE,scales = list(x = list(log = 10), y = list(log = 10)))
xyplot(pVal_ks ~ pVal_wilcox,ttRes,scales = list(x = list(log = 10), y = list(log = 10)))

# Binary matrix
ttResBin <- ttRes
ttResBin[ttRes < 0.5] <- TRUE
ttResBin[ttRes >= 0.5] <- FALSE
ttResBin

# How many pass threshold for each test?
testCounts <- colSums(ttResBin[,c("pVal_t","pVal_wilcox","pVal_ks")])
testCounts

# On how many tests does each gene pass the threshold?
geneCounts <- rowSums(ttResBin[,c("pVal_t","pVal_wilcox","pVal_ks")])
geneCounts

# How many genes pass by 3, 2, 1, 0 tests?
table(geneCounts)


