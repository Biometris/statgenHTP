
setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")

#_________________ EXAMPLE 1 _____________________________________________________________

#### DATASET-SPECIFIC FORMATING #######

inDat <- read.table("~/Documents/PostDoc/DROPS_Metabo/dataMetabo_ZB13_forSpatialModel_long_statgenHTP.csv",
                    h=T,sep=",")


#### BEGINING FUNCTION USE #######

inTP <- createTimePoints(dat = inDat[inDat$condition=="WW",],
                         genotype = "hybrid",
                         timePoint = "time",
                         repId = "Rep",
                         plotId = "pot",
                         rowNum = "y", colNum = "x",
                         addCheck = FALSE)

plot(inTP, plotType = "layout",
     timePoints = 1,
     highlight = c("B73_H","Oh43_H","PH207_H","B14a_H","Mo17_H"))

plot(inTP, plotType = "cor",
     traits = "value",
     timePoints = c(1:100))

plot(inTP, plotType = "box",
     traits = "value",
     timePoints = c(1:10))#,
     colorBy = "repId")

plot(inTP, plotType = "raw",
     traits = "value",
     timePoints = c(1:10), #c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     genotypes = c("B73_H","Oh43_H","PH207_H","B14a_H","Mo17_H"))

pdf("Phenovator_Rene_raw_data_na.pdf", height = 8, width = 12)
dev.off()


fitMods <- fitModels(TP = inTP,
                     trait = "value",
                     timePoints = c(1:10),
                     covariates = c("repId"))


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
plot(fitMods, plotType = "herit")
plot(fitMods, plotType = "variance")
plot(fitMods, plotType = "effDim")
plot(fitMods, plotType = "effDim")
plot(fitMods, plotType = "timeLapse", outFile = "spatialtrends.gif")

