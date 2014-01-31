# Care and feeding of data

gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderDataFiveYear.txt"
gDat <- read.delim(file = gdURL)
a <- list(veg = c("cabbage", "eggplant"), tNum = c(pi, exp(1), sqrt(2)), myAim = TRUE,joeNum = 2:6)
gdURL <- "http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderDataFiveYear.txt"
gDat <- read.delim(file = gdURL)
str(gdata)
str(gDat)
head(gDat)
tail(gDat)
peek(gDat)
names(gDat)
rownames(gDat)
head(rownames(gDat))
dim(gDat)
library(lattice)
xyplot(lifeExp ~ year, data=gDat)
xyplot(lifeExp ~ gdpPercap, data=gDat)
xyplot(lifeExp ~ gdpPercap, data=gDat,subset = country=="Colombia")
xyplot(lifeExp ~ gdpPercap, data=gDat,subset = country=="Colombia",type = c("p","r"))
xyplot(lifeExp ~ gdpPercap | continent, data=gDat,subset = country=="Colombia",type = c("p","r"))
xyplot(lifeExp ~ gdpPercap, data=gDat,group = continent, subset = year == 2007, auto.key=TRUE)
str(gDat)
head(gDat$lifeExp)
densityplot(~lifeExp,gDat)
table(gDat$lifeExp)
table(gDat$year)
class(gDat$year)
class(gDat$continent)
class(gDat$lifeExp)
summary(gDat$continent)
names(gDat$continent)
colnames(gDat$continent)
str(gDat$continent)
levels(gDat$continent)
barchart(gDat$continent)
dotplot(gDat$continent)
dotplot(gDat$continent,type="j")
dotplot(gDat$continent,type="h")
subset(gDat,subset=country=="Uruguay")
subset(gDat,subset=country=="Uruguay",select = c("country","year"))
with(subset(gDat,subset=country=="Uruguay",select = c("country","year")), mean(year))
with(subset(gDat,subset=country=="Uruguay"), mean(year))
     
# R objects and indexing
     
x <- 7
x[2] <- 100
x
x[1]
x[0]
y <- 1:3
z <- 3:7
y + z
n <- y + z
n
(x <- c("cabbage", pi, TRUE, 4.3))
a <- list("cabbage", pi, TRUE, 4.3)
a
names(a) <- c("veg", "dessert", "myAim", "number")
a <- list(veg = c("cabbage", "eggplant"), tNum = c(pi, exp(1), sqrt(2)), myAim = TRUE,
joeNum = 2:6))
a <- list(veg = c("cabbage", "eggplant"), tNum = c(pi, exp(1), sqrt(2)), myAim = TRUE,
joeNum = 2:6)
a
a[[2]]
a[2]
a$joeNum
a[c("tNum", "veg")]
n <- 8
set.seed(1)
(jDat <- data.frame(w = round(rnorm(n), 2),
x = 1:n,
y = I(LETTERS[1:n]),
z = runif(n) > 0.3,
v = rep(LETTERS[9:12], each = 2)))
w
jDat
jDat["w"]
jDat$w
jDat[["w"]]
w <- jDat$w
w < 0
which(w <0)
w[w<0]
w[seq(from = 1, to = length(w), by =2)]
w[-c(2,5)]
list("cabbage",pi,TRUE,4.3)
a <- list("cabbage",pi,TRUE,4.3)
a
names(a) <- c("veg", "dessert", "myAim","number")
a
a <- list(veg = c("cabbage", "eggplant"), tNum = c(pi, exp(1), sqrt(2)), myAim = TRUE,joeNum = 2:6)
a
a$joeNum
a[2]
a[[2]]
a["tNum"]
a[["tNum"]]
a[c("joeNum","veg")]
a[c("joeNum","veg")]
jMat <- outer(as.character(1:4), as.character(1:4), function(x, y) {
paste0("x", x, y)
})
jMat
rownames(jMat) <- paste0("row",seq_len(nrow(jMat)))
colnames(jMat) <- paste0("col",seq_len(nrow(jMat)))
dimnames(jMat)
jMat[2,3]
jMat[,3]
jMat[,3,drop=FALSE]
is.vector(jMat[,3])
is.vector(jMat[,3,drop=FALSE])
is.list(jMat[,3,drop=FALSE])
class(jMat[,3,drop=FALSE])
jMat[7]
jMat[c("row1","row2"),c("col1")]
jMat[c("row1","row2"),c("col1"),drop=FALSE]
jMat["row1",2:3] <- c("HEY!","THIS IS NUTS!")
jMat
matrix(1:15,nrow=5)
m<- matrix(c("yo!","foo?",nrow=3,ncol=6))
m
m<- matrix(c("yo!","foo?"),nrow=3,ncol=6,byrow=TRUE)
m

