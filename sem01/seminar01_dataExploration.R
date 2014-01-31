prDat <- read.table("GSE4051_MINI.txt",header=TRUE,row.names=1)
str(prDat)

# How many rows are there?
nrow(prDat)

# How many columns or variables are there?
ncol(prDat)

# Inspect the first few observations or the last few or a random sample
head(prDat)
tail(prDat)
prDat[sample(6),]

# What does row correspond to -- different genes or different mice?
# mice

# What are the variable names?
names(prDat)

# What "flavor" is each variable?
str(prDat)
# numeric, factor, factor, numeric, numeric, numeric

# For sample, make sure each integer between 1 and the number of rows occurs exactly once
samp <- sample(nrow(prDat))
x <- seq(1:nrow(prDat))
identical(sort(samp),x)

# For each factor variable, what are the levels?
levels(prDat$devStage)
levels(prDat$gType)

# How many observations do we have for each level of devStage? For gType?
table(prDat$devStage)
table(prDat$gType)

# Perform a cross-tabulation of devStage and gType
table(prDat$devStage,prDat$gType)

# If you had to take a wild guess, what do you think the intended experimental design was? What actually happened in real life?
# The effect of a knockout mutation on the expression of three genes

# For each quantitative variable, what are the extremes? How about average or median?
summary(prDat)

# Create a new data.frame called weeDat only containing observations for which expression of poisonFang is above 7.5
weeDat <- subset(prDat,subset= poisonFang > 7.5)

# For how many observations is poisonFang > 7.5? How do they break down by genotype and developmental stage?
nrow(weeDat)
table(weeDat$devStage,weeDat$gType)

# Print the observations with row names "Sample_16" and "Sample_38" to screen, showing only the 3 gene expression variables.
print(prDat[c("Sample_16","Sample_38"),c("crabHammer","eggBomb","poisonFang")])

# Which samples have expression of eggBomb less than the 0.10 quantile?
rownames(subset(prDat, subset= eggBomb < quantile(prDat$eggBomb,0.1)))