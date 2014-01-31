# ***** Take home - lattice *****

# Read in data and design file
library(lattice)
prDes <- read.table("GSE4051_design.tsv", header = TRUE, as.is = 1)
prDat <- read.table("GSE4051_data.tsv")

# Generate subsets of 2, 20, and 200
samp2 <- sample(1:nrow(prDat), size = 2, replace = FALSE)
hDat2 <- prDat[samp2,]

samp20 <- sample(1:nrow(prDat), size = 20, replace = FALSE)
hDat20 <- prDat[samp20,]

samp200 <- sample(1:nrow(prDat), size = 200, replace = FALSE)
hDat200 <- prDat[samp200,]

# Merge with design table for easier analysis
hDatDes2 <- merge(t(hDat2),prDes,by.x="row.names",by.y="sidChar")
hDatDes20 <- merge(t(hDat20),prDes,by.x="row.names",by.y="sidChar")
hDatDes200 <- merge(t(hDat200),prDes,by.x="row.names",by.y="sidChar")

# Reshape 2 gene data (made to be reusable)
oDat2 <-
  with(hDatDes2,
       data.frame(Row.names, sidNum, devStage, gType,
                  probeset = factor(rep(grep(".*at",colnames(hDatDes2),value=TRUE), each = nrow(hDatDes2))),
                  geneExp = c(unlist(hDatDes2[,grepl(".*at",colnames(hDatDes2))]))))

# Make a stripplot highlighting diffs between genes
stripplot(geneExp ~ devStage | probeset, oDat2, 
          layout = c(nlevels(oDat2$probeset), 1),
          groups = gType, auto.key = TRUE, jitter.data = TRUE)

# Make a stripplot highlighting diffs between groups
stripplot(geneExp ~ devStage | gType, oDat2, 
          layout = c(nlevels(oDat2$gType), 1),
          groups = probeset, auto.key = TRUE, jitter.data = TRUE)

# Do the same thing with density plots
densityplot(~ oDat2$geneExp, oDat2,
            groups = probeset, auto.key = TRUE)
densityplot(~ geneExp, oDat2,
            groups = gType, auto.key = TRUE)
bwplot(oDat2$geneExp ~ devStage | probeset, oDat2)
bwplot(oDat2$geneExp ~ devStage | gType, oDat2)

# Now let's make a heatmap with 20 gene data
library(RColorBrewer)
jPuBuGnFun <- colorRampPalette(brewer.pal(n=7,"PuBuGn"))
heatmap(as.matrix(hDat20),col=jPuBuGnFun(256))

# Too cute
pairs(hDat20[,sample(1:ncol(hDat20),size=6)])

# Eyesore!
splom(hDat20[,sample(1:ncol(hDat20),size=6)],panel=panel.smoothScatter,raster=TRUE)

# Just why...
library(hexbin)
hexplom(hDat20[,sample(1:ncol(hDat20),size=6)])