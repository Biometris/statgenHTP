
setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")

#_________________ EXAMPLE 1 _____________________________________________________________

#### DATASET-SPECIFIC FORMATING #######
source("./data-raw/exampleData.R")

# inDat <- data.table::fread("./data-raw/Original_PAM_reshape.csv",
                           # data.table = FALSE)

#### BEGINING FUNCTION USE #######

inTP <- createTimePoints(dat = PhenovatorDat1,
                         genotype = "Genotype",
                         timePoint = "timepoints",
                         repId = "Sowing_Block",
                         plotId = "pos",
                         rowNum = "y",
                         colNum = "x",
                         addCheck = TRUE,
                         checkGenotypes = c("col", "ely", "evo1", "ler"))

timepoint <- data.frame(attr(inTP,"timePoints"))
timepoint$timePoint <- lubridate::ymd_hms(timepoint$timePoint)

plot(inTP, plotType = "layout",
     timePoints = c(1),#c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     highlight = c("col", "ely", "evo1", "ler"))

plot(inTP, plotType = "cor",
     traits = "pheno",
     timePoints = c(1:10)) #c("2018-05-31 16:37:00", "2018-06-01 09:07:00"))

plot(inTP, plotType = "box",
     traits = "pheno",
     timePoints = c(1,3), #c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     colorBy = "repId")

plot(inTP, plotType = "raw",
     traits = "pheno",
     timePoints = c(1:10), #c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     genotypes = c("col", "ely", "evo1", "ler","1","18","340","2231"))


plot(inTP, plotType = "raw",
     traits = "pheno",
     # timePoints = c(1,25,60), #c("2018-05-31 16:37:00", "2018-06-01 09:07:00"),
     genotypes = c("col", "ely", "evo1", "ler","1","18","340","2231"))

plot(inTP, plotType = "layout",
     timePoints = c(1,3),
     highlight = c("col", "ely", "evo1", "ler"),
     outFile = "test.pdf")

plot(inTP, plotType = "cor",
     traits = "pheno",
     timePoints = c(1,4,14,23,42,50,68))

plot(inTP, plotType = "box", traits = "pheno", timePoints = c(5,6,7),
     colorBy = "repId")

plot(inTP, plotType = "raw", traits = "pheno",
     timePoints = c(8,4,2),
     genotypes = c("col", "ely", "evo1", "ler"),
     outFile = "test.pdf")

pdf("Phenovator_Rene_raw_data_na.pdf", height = 8, width = 12)
plot(inTP, plotType = "raw", traits = "pheno")
dev.off()

plot(inTP, plotType = "raw", traits = "pheno",
     genotypes = c("col", "ely", "evo1", "ler"))

fitMods <- fitModels(TP = inTP,
                     trait = "pheno",
                     covariates = c("repId", "Image_pos"),
                     useCheck = TRUE)

fitMods1b <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml")

fitMods1c <- fitModels(TP = inTP[1:3], trait = "pheno",
                       covariates = c("repId", "Image_pos"), engine = "asreml",
                       useCheck = TRUE)


spatCorr <- getCorrected(fitMods) #), outFile = "Corrected_PAM_modRep.csv")
spatCorr$timeNumber <- timepoint$timeNumber[match(spatCorr$timePoint,timepoint$timePoint)]

write.table(spatCorr,
            file="~/Documents/PostDoc/EPPN2020/Platform/Phenovator/Rene/spline_function_HTP/Corrected_PAM_modRepCheck.csv",
            row.names=F,sep=",")

genoPreds <- getGenoPred(fitMods, outFile = "BLUPs_PAM_modRep.csv")
variance <- getVar(fitMods)
h2 <- getHerit(fitMods)

genoPreds1b <- getGenoPred(fitMods1b)
genoPreds1c <- getGenoPred(fitMods1c)


plot(fitMods, plotType = "corrPred",
     genotypes = c("col", "ely", "evo1", "ler"), outFile = "test.pdf")
plot(fitMods, plotType = "rawPred")
plot(fitMods, plotType = "herit", timePoints = c(1,5))
plot(fitMods, plotType = "variance", outFile = "test.pdf")
plot(fitMods, plotType = "effDim")
plot(fitMods, plotType = "timeLapse", outFile = "spatialtrends.gif")

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
                          # repId = "Rep")

plot(inTP2, plotType = "layout",
     timePoints = c(4),
     highlight = c("GenoA1", "GenoA50", "GenoB1", "GenoB2") )

plot(inTP2, plotType = "cor",
     traits = "LA_Estimated")

plot(inTP2, plotType = "box",
     traits = "LA_Estimated",
     # timePoints = c(1:3),
     colorBy = "Scenario")

plot(inTP2, plotType = "raw",
     traits = "LA_Estimated",
     geno.decomp = "Scenario",
     timePoints = c(8,20,30),
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))

plot(inTP2, plotType = "raw",
     traits = "LA_Estimated",
     geno.decomp = "Scenario",
     timePoints = c(1:4),
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))


pdf("Phenovator_ZA17_raw_data_na.pdf", height = 8, width = 12)
plot(inTP2, plotType = "raw", traits = "LA_Estimated", geno.decomp = "Scenario")
dev.off()

fitMods2 <- fitModels(TP = inTP2[10:15],
                      trait = "LA_Estimated",
                      geno.decomp = c("Scenario", "population"))

is.null(fitMods2[[1]]$call$data$genotype)
is.null(fitMods2[[1]]$data$genotype)
is.null(fitMods2[[1]]$data)
fitMods2[[1]]

fitMods2b <- fitModels(TP = inTP2[10:15],
                       trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml")


fitMods2c <- fitModels(TP = inTP2[10:15],
                       trait = "LA_Estimated",
                       geno.decomp = c("Scenario", "population"),
                       engine = "asreml",
                       spatial = TRUE)

is.null(fitMods2c[[1]]$call$data$genotype)

is.null(fitMods2c[[1]]$data)

any(grep("SpATS",fitMods2[[1]]$call))
any(grep("SpATS",fitMods2c[[1]]$call))

fitMods2[[1]]$model$geno$geno.decomp
fitMods2c[[1]]$factor.names[1]
fitMods2[[1]]$model$geno$geno.decomp
fitMods2c[[1]]$factor.names[1]

## This crashes:
fitMods2crash <- fitModels(TP = inTP2[1:3],
                           trait = "LA_Estimated",
                           geno.decomp = "Scenario",
                           covariates = "population")

genoPreds2 <- getGenoPred(fitMods2, outFile = "BLUPs_ZA17_LeafArea.csv")
genoPreds2b <- getGenoPred(fitMods2b)
genoPreds2c <- getGenoPred(fitMods2c)

genoBLUPs2 <- getBLUPsGeno(fitMods2)
rowPreds2c <- getRowPred(fitMods2c)

spatCorr2 <- getCorrected(fitMods2, outFile = "Corrected_ZA17_LeafArea.csv")
variance2 <- getVar(fitMods2)
h22 <- getHerit(fitMods2)

#### doesn't work because the match is made for the geno while the objects use geno.decomp
#### work now !

plot(fitMods2,
     plotType = "corrPred",
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))

plot(fitMods2,
     plotType = "rawPred",
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))


## still not working because of other functions...
plot(fitMods2c,
     plotType = "rawPred",
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))

plot(fitMods2c,
     plotType = "corrPred",
     genotypes = c("GenoA4", "GenoA32", "GenoB32", "GenoB3"))


plot(fitMods2, plotType = "herit")

plot(fitMods2, plotType = "variance")
plot(fitMods2, plotType = "effDim")
plot(fitMods2, plotType = "timeLapse",
     outFile = "spatialtrends2.gif")

