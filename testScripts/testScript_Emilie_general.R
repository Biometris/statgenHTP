
setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")
library(statgenHTP)

#----------------------------------------------------------------------------------------
#####
##### Create a TP object containing the data from the Phenovator #######################
#####
data("PhenovatorDat1")
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Sowing_Block",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("col", "ely", "evo1", "ler"))

### Extract the time points table
timepoint <- data.frame(attr(phenoTP,"timePoints"))
timepoint$timePoint <- lubridate::ymd_hms(timepoint$timePoint)
### Summary of the TP object
# TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO

#----------------------------------------------------------------------------------------
#####
##### Plotting the TP object containing the data from the Phenovator #####################
#####

## ----layoutPlot
plot(phenoTP,
     plotType = "layout",
     timePoints = c(3))

## Plot the layout for the 3rd time point with the check genotypes highlighted.
plot(phenoTP,
     plotType = "layout",
     timePoints = c(3),
     highlight = c("col", "ely", "evo1", "ler"))

## Plot the layout for the 3rd time points.
plot(phenoTP,
     plotType = "layout",
     timePoints = c(3),
     highlight = c("col", "ely", "evo1", "ler"),
     showGeno = TRUE)

## Create a boxplot for PSII using the default all time points.
plot(phenoTP,
     plotType = "box",
     traits = "pheno")

## Create a boxplot for PSII with 5 time points and boxes colored by repIds within
## time point.
plot(phenoTP,
     plotType = "box",
     traits = "pheno",
     timePoints = c(1:5),
     colorBy = "repId")

## Create a boxplot for PSII with 5 time points and boxes grouped by repIds.
plot(phenoTP,
     plotType = "box",
     traits = "pheno",
     timePoints = c(1:5),
     groupBy = "repId")

## Create a correlation plot for PSII for a selection of time points.
plot(phenoTP,
     plotType = "cor",
     traits = "pheno",
     timePoints = seq(from=1,to=73,by=5))
#####

#----------------------------------------------------------------------------------------
#####-
##### Fitting SpATS models on the Phenovator data #######################################-
#####-

#####
##### ------ Fit a model few time points with geno = R ----------------------------------#
#####
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "pheno",
                        timePoints = seq(1,73,by=5)) #c(3,10,20,30,60)

obj <- modPhenoSp
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 6)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 6)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     # genotypes = c("col", "ler", "ely", "evo1","toto")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
     genotypes = c("toto", "tutu", "tata")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit", yLim = c(0,1))
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")


#####
##### ------ Fit a model few time points with geno = R and Covar ------------------------#
#####
modPhenoSpCov <- fitModels(TP = phenoTP,
                           trait = "pheno",
                           covariates = c("repId", "Image_pos"),
                           timePoints = seq(1,73,by=5))
obj <- modPhenoSpCov
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 6)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 6)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit")
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")

#####
##### ------ Fit a model few time points with geno = R and check genotypes --------------#
#####
modPhenoSpCheck <- fitModels(TP = phenoTP,
                             trait = "pheno",
                             useCheck = TRUE,
                             timePoints = seq(1,73,by=5))

modPhenoSpCheck <- fitModels(TP = phenoTP,
                             trait = "pheno",
                             useCheck = TRUE,
                             engine = "asreml",
                             spatial = TRUE,
                             timePoints = seq(1,73,by=5))
obj <- modPhenoSpCheck
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 6)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 6)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit", yLim = .5)
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")


#####
##### ------ Fit a model few time points with geno = R and rowcol design ----------------#
#####
modPhenoSpRCD <- fitModels(TP = phenoTP,
                            trait = "pheno",
                            timePoints = seq(1,73,by=20),
                            useRepId = TRUE)
obj <- modPhenoSpRCD
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 21)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 21)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit")
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")

#####
##### ------ Fit a model few time points with geno = R and rowcol design and cov Rep ----#
#####
modPhenoSpRCDCov <- fitModels(TP = phenoTP,
                              trait = "pheno",
                              covariates = c("repId","Image_pos"),
                              timePoints = seq(1,73,by=20),
                              useRepId = TRUE)
# Should display a warning? «repID already included as fixed effect when useRepId = TRUE»
summary(modPhenoSpRCDCov$`2018-05-31 16:37:00`)
obj <- modPhenoSpRCDCov

#####
##### ------ Fit a model few time points with geno = R, check = T, rowcol design and cov #
#####
modPhenoSpRCDCovCheck <- fitModels(TP = phenoTP,
                                   trait = "pheno",
                                   covariates = c("Image_pos"),
                                   useCheck = TRUE,
                                   timePoints = seq(1,73,by=20),
                                   useRepId = TRUE)

summary(modPhenoSpRCDCovCheck$`2018-06-05 16:37:00`)
obj <- modPhenoSpRCDCovCheck
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 21)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 21)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1"))
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit")
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")

#####
##### ------ Fit a model few time points with geno = F ----------------------------------#
#####
modPhenoSpFix <- fitModels(TP = phenoTP,
                           trait = "pheno",
                           timePoints = seq(1,73,by=10),
                           what="fixed")
##### geno = F, cov = image
modPhenoSpFix2 <- fitModels(TP = phenoTP,
                           trait = "pheno",
                           covariates = c("Image_pos"),
                           timePoints = seq(1,73,by=10),
                           what="fixed")
##### geno = F, cov = image, useRepId = T
modPhenoSpFix3 <- fitModels(TP = phenoTP,
                            trait = "pheno",
                            covariates = c("Image_pos"),
                            useRepId = TRUE,
                            timePoints = seq(1,73,by=10),
                            what="fixed")

obj <- modPhenoSpFix
obj <- modPhenoSpFix2
obj <- modPhenoSpFix3
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 11)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 11)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
effDim   <- getEffDims(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit")
## plot ED
plot(obj,
     plotType = "effDim",
     whichED = c("colId", "rowId", "fColRow","colfRow", "surface"))
## plot Var
plot(obj,
     plotType = "variance")


#----------------------------------------------------------------------------------------
#####-
##### Fitting ASReml models on the Phenovator data #######################################-
#####-






#####

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

