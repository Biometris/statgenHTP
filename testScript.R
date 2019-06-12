setwd(paste0("./output/"))

inDat <- data.table::fread("../data-raw/Original_PAM_reshape.csv",
                           data.table = FALSE)

# creating a unique ID per plant using the row and col coordinate
# (in principle the column "Sowing_Position" was also a unique ID but I wanted to see the position)
inDat$ID <- interaction(inDat[["x"]], inDat[["y"]], sep = "_")
inDat <- inDat[!is.na(inDat$pheno), ]
# Create an indicator for each plot (according to the row and column position)
inDat$pos <- paste0("c", inDat[["x"]], "r", inDat[["y"]])

# I removed a plant that has very few measurements
inDat <- inDat[inDat$pos != "c1r54",]

inTP <- createTimePoints(dat = inDat, genotype = "Genotype",
                         timePoint = "timepoints", repId = "Sowing_Block",
                         plotId = "pos", rowNum = "y", colNum = "x",
                         addCheck = TRUE,
                         checkGenotypes = c("col", "ely", "evo1", "ler"))

plot(inTP, plotType = "layout",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     highlight = c("col", "ely", "evo1", "ler"))
plot(inTP, plotType = "cor", traits = "pheno",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"))
plot(inTP, plotType = "box", traits = "pheno",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     colorBy = "repId")
plot(inTP, plotType = "raw", traits = "pheno",
     timePoints = c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     genotypes = c("col", "ely", "evo1", "ler"))

plot(inTP, plotType = "layout",
     timePoints = c(1,3),
     highlight = c("col", "ely", "evo1", "ler"), outFile = "test.pdf")
plot(inTP, plotType = "cor", traits = "pheno",  timePoints = c(4,23))
plot(inTP, plotType = "box", traits = "pheno", timePoints = c(5,6,7),
     colorBy = "repId")
plot(inTP, plotType = "raw", traits = "pheno",
     timePoints = c(8,4,2), genotypes = c("col", "ely", "evo1", "ler"),
     outFile = "test.pdf")

pdf("Phenovator_Rene_raw_data_na.pdf", height = 8, width = 12)
plot(inTP, plotType = "raw", traits = "pheno")
dev.off()

plot(inTP, plotType = "raw", traits = "pheno",
     genotypes = c("col", "ely", "evo1", "ler"))

fitMods <- fitModels(TP = inTP, trait = "pheno",
                     covariates = c("repId", "Image_pos"), useCheck = TRUE)
fitMods1b <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml")

fitMods1c <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       useCheck = TRUE)

genoPreds <- getGenoPred(fitMods, outFile = "BLUPs_PAM_modRep.csv")
colPreds <- getColPred(fitMods)
genoPreds1b <- getGenoPred(fitMods1b)
genoPreds1c <- getGenoPred(fitMods1c)
colPreds1c <- getColPred(fitMods1c)
spatCorr <- getCorrected(fitMods, outFile = "Corrected_PAM_modRep.csv")
variance <- getVar(fitMods)
h2 <- getHerit(fitMods)
genoBLUPS <- getBLUPsGeno(fitMods)

plot(fitMods, plotType = "corrPred",
     genotypes = c("col", "ely", "evo1", "ler"), outFile = "test.pdf")
plot(fitMods, plotType = "rawPred")
plot(fitMods, plotType = "herit", timePoints = c(1,5))
plot(fitMods, plotType = "variance", outFile = "test.pdf")
plot(fitMods, plotType = "effDim")
plot(fitMods, plotType = "timeLapse", outFile = "spatialtrends.gif")

## Second example
inDat2 <- data.table::fread("../data-raw/Data_modif_ZA17_anonymous.csv",
                            data.table = FALSE)
# Create an indicator for each plot (according to the row and column position)
inDat2$pos <- paste0("c", inDat2$Line, "r", inDat2$Position)

inTP2 <- createTimePoints(dat = inDat2, genotype = "geno", timePoint = "Date",
                          plotId = "pos", rowNum = "Position", colNum = "Line")

plot(inTP2, plotType = "layout", timePoints = "2017-04-13")
plot(inTP2, plotType = "cor", traits = "LA_Estimated")
plot(inTP2, plotType = "box", traits = "LA_Estimated",
     timePoints = names(inTP2[1:3]))
pdf("Phenovator_ZA17_raw_data_na.pdf", height = 8, width = 12)
plot(inTP2, plotType = "raw", traits = "LA_Estimated", geno.decomp = "Scenario")
dev.off()

fitMods2 <- fitModels(TP = inTP2[10:10], trait = "LA_Estimated",
                      geno.decomp = c("Scenario", "population"))
fitMods2b <- fitModels(TP = inTP2[10:10], trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml")
fitMods2c <- fitModels(TP = inTP2[10:10], trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml", spatial = TRUE)


## This crashes:
fitMods2crash <- fitModels(TP = inTP2[1:3], trait = "LA_Estimated",
                       geno.decomp = "Scenario", covariates = "population")

genoPreds2 <- getGenoPred(fitMods2, outFile = "BLUPs_ZA17_LeafArea.csv")
genoPreds2b <- getGenoPred(fitMods2b)
genoPreds2c <- getGenoPred(fitMods2c)

genoBLUPs2 <- getBLUPsGeno(fitMods2)
rowPreds2c <- getRowPred(fitMods2c)

spatCorr2 <- getCorrected(fitMods2, outFile = "Corrected_ZA17_LeafArea.csv")
variance2 <- getVar(fitMods2)
h22 <- getHerit(fitMods2)

plot(fitMods2, plotType = "corrPred")
plot(fitMods2, plotType = "rawPred")
plot(fitMods2, plotType = "herit")
plot(fitMods2, plotType = "variance")
plot(fitMods2, plotType = "effDim")
plot(fitMods2, plotType = "timeLapse", outFile = "spatialtrends2.gif")

