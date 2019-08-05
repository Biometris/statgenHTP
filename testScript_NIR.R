library(statgenHTP)
library(ggplot2)

inDat <- data.table::fread("./inst/extdata/Original_NIR_reshape.csv",
                           data.table = FALSE)
inDat$pos <- paste0("c", inDat[["x"]], "r", inDat[["y"]])

inTP <- createTimePoints(inDat, genotype = "Genotype", timePoint = "timepoints",
                         plotId = "pos", repId = "Replicate", rowNum = "y",
                         colNum = "x", addCheck = TRUE,
                         checkGenotypes = c("col", "ely", "evo1", "ler"))

plot(inTP, timePoints = 1)
plot(inTP, plotType = "box", traits = "pheno")
plot(inTP, plotType = "cor", traits = "pheno")
plot(inTP, plotType = "raw", traits = "pheno")
plot(inTP, plotType = "raw", traits = "pheno",
     genotypes = c("col", "ely", "evo1", "ler"))

fitMods1a <- fitModels(TP = inTP, trait = "pheno",
                       covariates = c("repId", "Image.pos"), useCheck = TRUE)

preds <- getGenoPred(fitMods1a)
corr <- getCorrected(fitMods1a)
effDim <- getEffDims(fitMods1a)
effRat <- getEffDims(fitMods1a, EDType = "ratio")
variance <- getVar(fitMods1a)
herit <- getHerit(fitMods1a)

plot(fitMods1a, plotType = "spatial", timePoints = 1)
plot(fitMods1a, plotType = "spatial", timePoints = 1, spaTrend = "percentage")
plot(fitMods1a, plotType = "timeLapse", outFile = "test_NIR.gif")
plot(fitMods1a, plotType = "rawPred", genotypes = c("col", "ely", "evo1", "ler"))
plot(fitMods1a, plotType = "corrPred", genotypes = c("col", "ely", "evo1", "ler"))
plot(fitMods1a, plotType = "herit")
plot(fitMods1a, plotType = "effDim")
plot(fitMods1a, plotType = "effDim", EDType = "ratio")
plot(fitMods1a, plotType = "variance")

fitMods1b <- fitModels(TP = inTP, trait = "pheno",
                       covariates = c("repId", "Image.pos"), useCheck = TRUE,
                       engine = "asreml", spatial = TRUE)
attr(fitMods1b[[1]], "sumTab")

preds2 <- getGenoPred(fitMods1b)
corr2 <- getCorrected(fitMods1b)
variance2 <- getVar(fitMods1b)
herit2 <- getHerit(fitMods1b)

plot(fitMods1b, plotType = "spatial", timePoints = 1)
plot(fitMods1b, plotType = "rawPred", genotypes = c("col", "ely", "evo1", "ler"))
plot(fitMods1b, plotType = "corrPred", genotypes = c("col", "ely", "evo1", "ler"))
plot(fitMods1b, plotType = "herit")
plot(fitMods1b, plotType = "variance")

compPred <- merge(preds, preds2, by = c("genotype", "timePoint"))
cor(compPred$predicted.values.x, compPred$predicted.values.y)
compPred$timePoint <- as.factor(compPred$timePoint)
ggplot(compPred, aes(x = predicted.values.x, y = predicted.values.y)) +
  geom_point() + geom_abline(colour = "red") +
  ggforce::facet_wrap_paginate(~timePoint, nrow = 5, ncol = 5, scales = "free")

compCorr <- merge(corr, corr2, by = c("plotId", "timePoint"))
cor(compCorr$newTrait.x, compCorr$newTrait.y)
compCorr$timePoint <- as.factor(compCorr$timePoint)
ggplot(compCorr, aes(x = newTrait.x, y = newTrait.y)) +
  geom_point() + geom_abline(colour = "red") +
  ggforce::facet_wrap_paginate(~timePoint, nrow = 5, ncol = 5, scales = "free")

compHerit <- merge(herit, herit2, by = "timePoint")
cor(compHerit$h2.x, compHerit$h2.y)
plot(compHerit$h2.x, compHerit$h2.y)
abline(0, 1)
