## All test cases Phenovator

library(statgenHTP)

# Create inTP ----
{
  inDat <- data.table::fread("./inst/extdata//Original_PAM_reshape.csv",
                             data.table = FALSE)
  inDat <- inDat[!is.na(inDat$pheno), ]
  # Create an indicator for each plot (according to the row and column position)
  inDat$pos <- paste0("c", inDat[["x"]], "r", inDat[["y"]])
  # I removed a plant that has very few measurements
  inDat <- inDat[inDat$pos != "c1r54",]
  inTP <- createTimePoints(dat = inDat, genotype = "Genotype",
                           timePoint = "timepoints",
                           repId = "Sowing_Block",
                           plotId = "pos",
                           rowNum = "y", colNum = "x",
                           addCheck = TRUE,
                           checkGenotypes = c("col", "ely", "evo1", "ler"))
}

# Define function for comparing outputs ----
compare <- function(fitMods1, fitMods2, outFile = "comp") {
  ## Compare spatial corrected values.
  spatCor1 <- getCorrected(fitMods1)
  spatCor2 <- getCorrected(fitMods2)
  compSp <- merge(spatCor1, spatCor2, by = c("plotId", "timePoint"))
  print(cor(compSp$newTrait.x, compSp$newTrait.y, use = "pair"))
  compSp$timePoint <- as.factor(compSp$timePoint)
  pdf(file.path("comparisons", paste0(outFile, "Sp.pdf")))
  for (i in seq_along(fitMods1)) {
    p <- ggplot2::ggplot(compSp, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
      ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
      ggforce::facet_wrap_paginate('timePoint', nrow = 1, ncol = 1,
                                   scales = "free", page = i)
    plot(p)
  }
  dev.off()
  ## Compare genotype predictions.
  genoPred1 <- getGenoPred(fitMods1)
  genoPred2 <- getGenoPred(fitMods2)
  compGP <- merge(genoPred1, genoPred2, by = c("genotype", "timePoint"))
  compGP$timePoint <- as.factor(compGP$timePoint)
  print(cor(compGP$predicted.values.x, compGP$predicted.values.y))
  pdf(file.path("comparisons", paste0(outFile, "GP.pdf")))
  for (i in seq_along(fitMods1)) {
    p <- ggplot2::ggplot(compGP, ggplot2::aes(x = predicted.values.x,
                                              y = predicted.values.y)) +
      ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
      ggforce::facet_wrap_paginate('timePoint', nrow = 1, ncol = 1,
                                   scales = "free", page = i)
    plot(p)
  }
  dev.off()
}

# Set time points ----
tp <- c(1, 11, 21, 31, 41, 51, 61, 71)

# From the original script ----
fitMods1a <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                       covariates = c("repId", "Image_pos"),
                       useCheck = TRUE)
fitMods1b <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                       covariates = c("repId", "Image_pos"),
                       useCheck = TRUE,
                       engine = "asreml", spatial = TRUE)
compare(fitMods1a, fitMods1b, outFile = "phenovatorOrig")

# Remove covariates ----
fitMods1_1a <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         useCheck = TRUE)
fitMods1_1b <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         useCheck = TRUE,
                         engine = "asreml", spatial = TRUE)
compare(fitMods1_1a, fitMods1_1b, outFile = "phenovatorNoCovar")

# Remove check ----
fitMods1_2a <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         covariates = c("repId", "Image_pos"))
fitMods1_2b <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         covariates = c("repId", "Image_pos"),
                         engine = "asreml", spatial = TRUE)
compare(fitMods1_2a, fitMods1_2b, outFile = "phenovatorNoCheck")

# Not spatial ----
fitMods1_3b <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         covariates = c("repId", "Image_pos"),
                         useCheck = TRUE,
                         engine = "asreml")
compare(fitMods1a, fitMods1_3b, outFile = "phenovatorNoSpatial")

compare(fitMods1a, fitMods1_1a, outFile = "phenovatorSpATS_withWithoutCovar")
compare(fitMods1b, fitMods1_1b, outFile = "phenovatorAsreml_withWithoutCovar")

# Add geno.decomp ----
fitMods1_4a <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         geno.decomp = "Basin")
fitMods1_4b <- fitModels(TP = inTP, trait = "pheno", timePoints = tp,
                         geno.decomp = "Basin",
                         engine = "asreml", spatial = TRUE)
compare(fitMods1_4a, fitMods1_4b, outFile = "phenovatorGenoDecomp")

# Compare asreml and lme4
pdf("./comparisons/phenovatorAsremlLme4.pdf")
for (i in tp) {
  testAs <- fitModels(TP = inTP, trait = "pheno", timePoints = i,
                      covariates = c("repId", "Image_pos"),
                      engine = "asreml")
  predAs <- getGenoPred(testAs)
  testLm <- lme4::lmer(pheno~repId+Image_pos+(1|genotype), inTP[[i]])
  predLm <- coef(testLm)$genotype[, 1, drop = FALSE] +
    mean(c(0, unlist(coef(testLm)$genotype[1,2:8]))) +
    mean(c(0, unlist(coef(testLm)$genotype[1,9:19])))
  plot(predLm$`(Intercept)`, predAs$predicted.values)
  abline(0, 1)
}
dev.off()
