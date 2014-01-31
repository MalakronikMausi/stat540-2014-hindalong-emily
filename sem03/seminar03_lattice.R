library(lattice)

# Read in mini data
kDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)

# Do stuff
table(kDat$gType)
with(kDat,table(devStage,gType))
xyplot(eggBomb ~ crabHammer, kDat)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       auto.key = TRUE)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       outer=TRUE)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       outer = TRUE, grid = TRUE)
xyplot(eggBomb + poisonFang ~ crabHammer, kDat,
       outer = TRUE, grid = TRUE,
       groups = gType, auto.key = TRUE)

nDat <-
  with(kDat,
       data.frame(sample, devStage, gType, crabHammer,
                  probeset = factor(rep(c("eggBomb", "poisonFang"), each = nrow(kDat))),
                  geneExp = c(eggBomb, poisonFang)))
str(nDat)

xyplot(geneExp ~ crabHammer | probeset, nDat,
       grid = TRUE,
       groups = gType, auto.key = TRUE)

oDat <-
  with(kDat,
       data.frame(sample, devStage, gType,
                  probeset = factor(rep(c("crabHammer", "eggBomb",
                                          "poisonFang"), each = nrow(kDat))),
                  geneExp = c(crabHammer, eggBomb, poisonFang)))
str(oDat)

stripplot(~ geneExp, oDat)
stripplot(~ geneExp | probeset, oDat)
stripplot(probeset ~ geneExp, oDat)
stripplot(probeset ~ geneExp | probeset, oDat)
stripplot(probeset ~ geneExp, oDat, jitter.dat = TRUE)

stripplot(~ geneExp | probeset, oDat,
          layout = c(nlevels(oDat$probeset), 1))

stripplot(~ geneExp| probeset, oDat, 
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE)

stripplot(~ geneExp| probeset, oDat, 
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE, jitter.data = TRUE)

stripplot(geneExp ~ devStage, oDat)

stripplot(geneExp ~ devStage | probeset, oDat, 
          layout = c(nlevels(oDat$probeset), 1))

stripplot(geneExp ~ devStage | probeset, oDat, 
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE)

stripplot(geneExp ~ devStage | probeset, oDat, 
          layout = c(nlevels(oDat$probeset), 1),
          groups = gType, auto.key = TRUE, grid = TRUE,
          type = c('p', 'a'))

densityplot(~ geneExp, oDat)

densityplot(~ geneExp | gType, oDat,
            grid = TRUE,layout = c(nlevels(oDat$gType), 1))

densityplot(~ geneExp, oDat,
            groups = gType, auto.key = TRUE)

jBw <- 0.2
jn <- 400
densityplot(~ geneExp, oDat,
            groups = gType, auto.key = TRUE,
            bw = jBw, n = jn,
            main = paste("bw =", jBw, ", n =", jn))

bwplot(geneExp ~ devStage, oDat)

bwplot(geneExp ~ devStage | gType, oDat)

bwplot(geneExp ~ devStage, oDat,
       panel = panel.violin)

# Read in big data and design file
prDes <- read.table("GSE4051_design.tsv", header = TRUE, as.is = 1)
prDat <- read.table("GSE4051_data.tsv")

yo <- sample(1:nrow(prDat), size = 50, replace = TRUE)
hDat <- prDat[yo, ]
str(hDat)

hDat <- as.matrix(t(hDat))
rownames(hDat) <- with(prDes,
                       paste(devStage, gType, sidChar, sep="_"))
str(hDat)

heatmap(hDat, Rowv = NA, Colv = NA, col = cm.colors(256),
        scale="none", margins = c(5, 8))

library(RColorBrewer)

jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))

heatmap(hDat, Rowv = NA, Colv = NA, scale="none", margins = c(5, 8),
        col = jBuPuFun(256))

heatmap(hDat, margins = c(5, 8), col = jBuPuFun(256))

heatmap(hDat, col = jBuPuFun(256), margins = c(5, 8), scale=c("column"))

library(gplots)
heatmap.2(hDat, col = jBuPuFun, trace = "none")

set.seed(924)
(yo <- sample(1:ncol(prDat), size = 2))
y <- prDat[[yo[1]]]
z <- prDat[[yo[2]]]

xyplot(y ~ z, asp=1)
xyplot(y ~ z, asp=2)

smoothScatter(y ~ z, asp = 1)

xyplot(y ~ z, asp = 1, panel = panel.smoothScatter, nbin = 10)
xyplot(y ~ z, asp = 1, panel = panel.smoothScatter, nbin = 150)

library(hexbin)
hexbinplot(y~z)

set.seed(924)
(yo <- sample(1:ncol(prDat), size = 4))
pairDat<-subset(prDat,select=yo)

pairs(pairDat)
pairs(pairDat,panel=function(...) smoothScatter(...,add=TRUE))
splom(pairDat)
splom(pairDat,panel=panel.smoothScatter,raster=TRUE)
hexplom(pairDat)