# Warm-up

rnorm(n = 10)
rnorm(n,mean=0,sd=1)
set.seed(540)
rnorm(10)

n <- 10
B <- 4
x <- matrix(rnorm(n * B), nrow = n)

rownames(x) <- sprintf("obs%02d",1:n)
colnames(x) <- sprintf("samp%02d",1:B)
dim(x)

mean(x[,2])
mean(apply(x,2,mean))

# Distributuion of sample means - normal distribution

numSamps = 100
a = matrix(rnorm(1 * numSamps), nrow = 1)
b = matrix(rnorm(10 * numSamps), nrow = 10)
c = matrix(rnorm(100 * numSamps), nrow = 100)
d = matrix(rnorm(1000 * numSamps), nrow = 1000)
e = matrix(rnorm(10000 * numSamps), nrow = 10000)
allSamples = list(a,b,c,d,e)

trueMeanSD <- c(1/sqrt(1),1/sqrt(10),1/sqrt(100),1/sqrt(1000),1/sqrt(10000))
sampMeanSD <- lapply(allSamples,function(x) {sd(colMeans(x))})
sampMeanIQR <- lapply(allSamples,function(x) {IQR(colMeans(x))})
sampMeanMad <- lapply(allSamples,function(x) {mad(colMeans(x))})

stats <- cbind(trueMeanSD,sampMeanSD,sampMeanIQR,sampMeanMad)
rownames(stats) <- c("n1","n10","n100","n1000","n10000")
colnames(stats) <- c("trueMeanSD","sampMeanSD","sampMeanIQR","sampMeanMed")

# Distributuion of sample means - binomial distribution

numSamps = 100
a = rbinom(n=numSamps, size=1, prob = 0.7) / 1
b = rbinom(n=numSamps, size=10, prob = 0.7) / 10
c = rbinom(n=numSamps, size=100, prob = 0.7) / 100
d = rbinom(n=numSamps, size=1000, prob = 0.7) / 1000
e = rbinom(n=numSamps, size=10000, prob = 0.7) / 10000
allSamples = list(a,b,c,d,e)

trueMeanSD <- c(0.3/sqrt(1),0.3/sqrt(10),0.3/sqrt(100),0.3/sqrt(1000),0.3/sqrt(10000))
sampMeanSD <- lapply(allSamples,sd)
sampMeanIQR <- lapply(allSamples,IQR)
sampMeanMad <- lapply(allSamples,mad)

stats <- cbind(trueMeanSD,sampMeanSD,sampMeanIQR,sampMeanMad)
rownames(stats) <- c("n1","n10","n100","n1000","n10000")
colnames(stats) <- c("trueMeanSD","sampMeanSD","sampMeanIQR","sampMeanMed")

# CDF

x <- rnorm(10,mean=100,sd=12)
difobserved <- length(x[x <= 90]) / 10
expected <- pnorm(90,mean=100,sd=12)
diff1 <- expected - observed

y <- rnorm(100,mean=100,sd=12)
observed <- length(y[y <= 90]) / 100
expected <- pnorm(90,mean=100,sd=12)
diff2 <- expected - observed

z <- rnorm(1000,mean=100,sd=12)
observed <- length(z[z <= 90]) / 1000
expected <- pnorm(90,mean=100,sd=12)
diff3 <- expected - observed

xb <- rbinom(n=100,size=10,prob=0.2)
observed <- length(xb[xb<=3]) / 100
expected <- pbinom(3,size=10,prob=0.2)
diffxb <- expected - observed

yb <- rbinom(n=100,size=100,prob=0.2)
observed <- length(xb[xb<=30]) / 100
expected <- pbinom(30,size=100,prob=0.2)
diffyb <- expected - observed

zb <- rbinom(n=100,size=1000,prob=0.2)
observed <- length(xb[xb<=300]) / 100
expected <- pbinom(300,size=1000,prob=0.2)
diffzb <- expected - observed

# Plotting

densityplot(x)
densityplot(y)
densityplot(z)
densityplot(xb)
densityplot(yb)
densityplot(zb)

# More available on website...