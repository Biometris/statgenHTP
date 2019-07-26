## All test cases phenoarch

library(statgenHTP)

# Create inTP2 ----
{
  ## Simple example.
  inDat2 <- data.table::fread("./inst/extdata/Data_modif_ZA17_anonymous.csv",
                              data.table = FALSE)
  # Create an indicator for each plot (according to the row and column position)
  inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)
  inTP2 <- createTimePoints(dat = inDat2,
                            genotype = "geno",
                            timePoint = "Date",
                            plotId = "pos",
                            rowNum = "Position",
                            colNum = "Line",
                            addCheck = TRUE,
                            checkGenotypes = c("GenoA1", "GenoB2", "GenoA3"))
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
tp <- c(1, 6, 11, 16, 21, 16, 31, 35)

# From the original script ----
fitMods2a <- fitModels(TP = inTP2, trait = "LA_Estimated", timePoints = tp,
                       geno.decomp = c("Scenario", "population"))
fitMods2b <- fitModels(TP = inTP2, trait = "LA_Estimated", timePoints = tp,
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml", spatial = TRUE)
compare(fitMods2a, fitMods2b, outFile = "phenoarchOrig")

# Add fake check ----
fitMods2_1a <- fitModels(TP = inTP2, trait = "LA_Estimated", timePoints = tp,
                         geno.decomp = c("Scenario", "population"),
                         useCheck = TRUE)
fitMods2_1b <- fitModels(TP = inTP2, trait = "LA_Estimated", timePoints = tp,
                         geno.decomp = c("Scenario", "population"),
                         useCheck = TRUE,
                         engine = "asreml", spatial = TRUE)
compare(fitMods2_1a, fitMods2_1b, outFile = "phenoarchCheck")



