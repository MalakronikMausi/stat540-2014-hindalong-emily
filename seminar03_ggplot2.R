library(ggplot2)

apropos("^geom_")
apropos("^stat_")
apropos("^aes_")
apropos("^scale_")

qplot(crabHammer, eggBomb, data = kDat)

p <- ggplot(kDat, aes(x = crabHammer, y = eggBomb))

(p  <- p + geom_point())

(p <- p + stat_smooth())

(p <- p + theme_bw() + 
   xlab("Expression of crabHammer") + 
   ylab("Expression of eggBomb") + 
   ggtitle("Scatterplot for expression levels"))



nDat <-
  with(kDat,
       data.frame(sample, devStage, gType, crabHammer,
                  probeset = factor(rep(c("eggBomb", "poisonFang"),
                                        each = nrow(kDat))),
                  geneExp = c(eggBomb, poisonFang)))
str(nDat)

(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point())

(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point() + 
   stat_smooth(se = F))

(p <- ggplot(nDat, aes(crabHammer, geneExp, color = probeset)) + 
   geom_point() + 
   stat_smooth(se = F, aes(group = 1)))

(p <- ggplot(nDat, aes(crabHammer, geneExp,color = probeset)) + 
   geom_point() + 
   facet_wrap(~ probeset))

(p <- ggplot(nDat, aes(crabHammer, geneExp, color = gType)) + 
   geom_point() + 
   facet_wrap(~ probeset))

# To be continued...