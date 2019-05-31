folder <- "C:/Projects/R packages/statgenhtp/"
setwd(paste0(folder, "output/"))

inDat <- data.table::fread("../data-raw/Original_PAM_reshape.csv",
                               data.table = FALSE)

# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat$x, "r", inDat$y)

# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]
inDat <- droplevels(inDat)

inTD <- createTD(dat = inDat, genotype = "Genotype",
                 timePoint = "timepoints", plotId = "pos", rowNum = "y",
                 colNum = "x", addCheck = TRUE,
                 checkGenotypes = c("col", "ely", "evo1", "ler"))

plot(inTD, plotType = "layout",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     highlight = c("col", "ely", "evo1", "ler"))
plot(inTD, plotType = "cor", traits = "pheno")
plot(inTD, plotType = "box", traits = "pheno",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     colorBy = "Sowing_Block")

pdf("Phenovator_Rene_raw_data_na.pdf", height = 8, width = 12)
plot(inTD, plotType = "raw", traits = "pheno")
dev.off()

fitMods <- fitModels(TD = inTD, trait = "pheno",
                     covariates = c("Sowing_Block", "Image_pos"))

genoPreds <- getGenoPred(fitMods, outFile = "BLUPs_PAM_modRep.csv")
spatCorr <- getCorrected(fitMods, outFile = "Corrected_PAM_modRep.csv")
variance <- getVar(fitMods)
h2 <- getHerit(fitMods)

plot(fitMods, plotType = "corrPred")
plot(fitMods, plotType = "rawPred")
plot(fitMods, plotType = "herit")
plot(fitMods, plotType = "variance")
plot(fitMods, plotType = "effDim")

## Second example
inDat2 <- data.table::fread("../data-raw/Data_modif_ZA17_anonymous.csv",
                           data.table = FALSE)

# Create an indicator for each plot (according to the row and column position)
inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)

inTD2 <- createTD(dat = inDat2, genotype = "geno", timePoint = "Date",
                  plotId = "pos", rowNum = "Position", colNum = "Line")

plot(inTD2, plotType = "layout", timePoints = "2017-04-13")
plot(inTD2, plotType = "cor", traits = "LA_Estimated")
plot(inTD2, plotType = "box", traits = "LA_Estimated",
     timePoints = names(inTD2[1:3]))
pdf("Phenovator_ZA17_raw_data_na.pdf", height = 8, width = 12)
plot(inTD2, plotType = "raw", traits = "LA_Estimated")
dev.off()

fitMods2 <- fitModels(TD = inTD2, trait = "LA_Estimated",
                      geno.decomp = "Scenario", covariates = "population")

genoPreds2 <- getGenoPred(fitMods2, outFile = "BLUPs_ZA17_LeafArea.csv")
spatCorr2 <- getCorrected(fitMods2, outFile = "Corrected_ZA17_LeafArea.csv")
variance2 <- getVar(fitMods2)
h22 <- getHerit(fitMods2)

plot(fitMods2, plotType = "corrPred")
plot(fitMods2, plotType = "rawPred")
plot(fitMods2, plotType = "herit")
plot(fitMods2, plotType = "variance")
plot(fitMods2, plotType = "effDim")


