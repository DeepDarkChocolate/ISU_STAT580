install.packages("teigen")
install.packages("clusterGeneration")
library(teigen)
library(clusterGeneration)


sim <- genRandomClust(2, sepVal = .35, numReplicate = 1,
                      outputDatFlag = FALSE, outputLogFlag = FALSE, outputEmpirical = FALSE,
                      outputInfo = FALSE)$datList[[1]]

par(mfrow = c(2, 2))

CIIC <- teigen(sim, 2, models = "CIIC", verbose = FALSE,
               dfupdate = "numeric")
plot(CIIC, levels = seq(0.03, 1, by = 0.1), what = "contour",
     xlab = "Variable 1", ylab = "Variable 2", main = "CIIC model")

irisknown <- iris[,5]
irisknown[134:150] <- NA
triris <- teigen(iris[,-5], models="CUUU", init="uniform", known=irisknown)
teigen
View(teigen)
summary(CIIC)
CIIC$iter
CIIC$parameters
sum(CIIC$iclresults$classification != CIIC$classification)
kmeans1 <- kmeans(sim, CIIC$G)
sum(kmeans1$cluster != CIIC$classification)
kmeans1

