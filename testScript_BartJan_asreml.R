
#setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")

#_________________ EXAMPLE 1 _____________________________________________________________

#### DATASET-SPECIFIC FORMATING #######

inDat <- data.table::fread("./data-raw/Original_PAM_reshape.csv",
                           data.table = FALSE)

# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat[["x"]], "r", inDat[["y"]])

# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]

#### BEGINING FUNCTION USE #######

inTP <- createTimePoints(dat = inDat, genotype = "Genotype",
                         timePoint = "timepoints",
                         repId = "Sowing_Block",
                         plotId = "pos",
                         rowNum = "y", colNum = "x",
                         addCheck = TRUE,
                         checkGenotypes = c("col", "ely", "evo1", "ler"))

fitMods <- fitModels(TP = inTP, trait = "pheno", timePoints = 1:3,
                     covariates = c("repId", "Image_pos"))
fitMods1a <- fitModels(TP = inTP, trait = "pheno", timePoints = 1:3,
                       covariates = c("repId", "Image_pos"), useCheck = TRUE)
fitMods1b <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml")
fitMods1c <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       useCheck = TRUE)
fitMods1d <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       spatial = TRUE)

fitMods1e <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       useCheck = TRUE, spatial = TRUE)

spatCorr <- getCorrected(fitMods)
spatCorr1a <- getCorrected(fitMods1a)
spatCorr1b <- getCorrected(fitMods1b)
spatCorr1c <- getCorrected(fitMods1c)
spatCorr1d <- getCorrected(fitMods1d)
spatCorr1e <- getCorrected(fitMods1e)

comp <- merge(spatCorr, spatCorr1b, by = c("timePoint", "plotId"))
comp$timePoint <- as.factor(comp$timePoint)
ggplot2::ggplot(comp, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
  ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
  ggplot2::facet_wrap(~timePoint)

comp1a <- merge(spatCorr1a, spatCorr1c, by = c("timePoint", "plotId"))
comp1a$timePoint <- as.factor(comp1a$timePoint)
ggplot2::ggplot(comp1a, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
  ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
  ggplot2::facet_wrap(~timePoint)

comp1b <- merge(spatCorr, spatCorr1d, by = c("timePoint", "plotId"))
comp1b$timePoint <- as.factor(comp1b$timePoint)
ggplot2::ggplot(comp1b, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
  ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
  ggplot2::facet_wrap(~timePoint)

comp1c <- merge(spatCorr1a, spatCorr1e, by = c("timePoint", "plotId"))
comp1c$timePoint <- as.factor(comp1c$timePoint)
ggplot2::ggplot(comp1c, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
  ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
  ggplot2::facet_wrap(~timePoint)

#_________________ EXAMPLE 2 _____________________________________________________________

#### DATASET-SPECIFIC FORMATING #######

inDat2 <- data.table::fread("./data-raw/Data_modif_ZA17_anonymous.csv",
                            data.table = FALSE)
# Create an indicator for each plot (according to the row and column position)
inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)

#### BEGINING FUNCTION USE #######

inTP2 <- createTimePoints(dat = inDat2,
                          genotype = "geno",
                          timePoint = "Date",
                          plotId = "pos",
                          rowNum = "Position",
                          colNum = "Line") #,


fitMods2 <- fitModels(TP = inTP2[13:15],
                      trait = "LA_Estimated",
                      geno.decomp = c("Scenario", "population"))

fitMods2b <- fitModels(TP = inTP2[13:15],
                       trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml")


fitMods2c <- fitModels(TP = inTP2[13:15],
                       trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml",
                       spatial = TRUE)

## This crashes:
fitMods2crash <- fitModels(TP = inTP2[1:3],
                           trait = "LA_Estimated",
                           geno.decomp = "Scenario",
                           covariates = "population")

spatCorr2 <- getCorrected(fitMods2)

## This crashes since there is nothing to correct for.
#spatCorr2b <- getCorrected(fitMods2b)

spatCorr2c <- getCorrected(fitMods2c)

# comp2 <- merge(spatCorr2, spatCorr2b, by = c("timePoint", "plotId"))
# comp2$timePoint <- as.factor(comp2$timePoint)
# ggplot2::ggplot(comp2, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
#   ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
#   ggplot2::facet_wrap(~timePoint)

comp2a <- merge(spatCorr2, spatCorr2c, by = c("timePoint", "plotId"))
comp2a$timePoint <- as.factor(comp2a$timePoint)
ggplot2::ggplot(comp2a, ggplot2::aes(x = newTrait.x, y = newTrait.y)) +
  ggplot2::geom_point() + ggplot2::geom_abline(colour = "red") +
  ggplot2::facet_wrap(~timePoint)



