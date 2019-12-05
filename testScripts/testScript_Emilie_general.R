
setwd(paste0("~/Documents/PostDoc/Pipeline_Biometris/statgenHTP/"))
# devtools::load_all(".")
library(statgenHTP)

#----------------------------------------------------------------------------------------
#####
##### Create a TP object containing the data from the Phenovator #######################
#####
data("PhenovatorDat1")
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Replicate",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("check1", "check2", "check3", "check4"))

### Extract the time points table
timepoint <- attr(phenoTP,"timePoints")
timepoint$timePoint <- lubridate::ymd_hms(timepoint$timePoint)
### Summary of the TP object
# TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO TO DO
scales::pretty_breaks(n=2,
                      min.n=1,
                      shrink.sml = 10)(timepoint$timeNumber)

floor(max(timepoint$timeNumber)/3)*c(1,2,3)
max(timepoint$timePoint)/3*c(1,2,3)

int <- lubridate::interval(timepoint$timePoint[1], timepoint$timePoint[max(timepoint$timeNumber)])
sec <- lubridate::int_length(int)
sec/4
timepoint$timePoint[1]+(sec/4*c(1,2,3))


#----------------------------------------------------------------------------------------
#####
##### Create a TP object containing the data from the Phenoarch #######################
#####
data("PhenoarchDat1")
phenoTParch <- createTimePoints(dat = PhenoarchDat1,
                                experimentName = "Phenoarch",
                                genotype = "geno",
                                timePoint = "Date",
                                repId = "Rep",
                                plotId = "pos",
                                rowNum = "Position",
                                colNum = "Line")

### Extract the time points table
timepointArch <- data.frame(attr(phenoTParch,"timePoints"))
timepointArch$timePoint <- lubridate::ymd(timepointArch$timePoint)
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
                        trait = "EffpsII",
                        timePoints = seq(from = 1, to = 73, by = 5))

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
# plot Spat
plot(modPhenoSp,
     timePoints = 36,
     plotType = "spatial",
     spaTrend = "percentage")
plot(modPhenoSp,
     timePoints = 36,
     plotType = "spatial",
     spaTrend = "raw")

# plot(modPhenoSp,
#      plotType = "spatial",
#      spaTrend = "percentage",
#      outFile = "spatial_percentage.pdf",
#      outFileOpts = c(width = 8.5, height = 5))

#####
##### ------ splines and outliers ------------------------#
#####
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        useCheck = TRUE)

spatCorr <- getCorrected(modPhenoSp)
write.table(spatCorr,
            "~/Documents/PostDoc/EPPN2020/Outliers/Doc_Git/time-course/Corrected_EffPsII.csv",
            sep=",", row.names=F)

fit.spline <- fitSpline(corrDat = spatCorr,
                        trait = "EffpsII",
                        knots = 50)

pred.Dat <- fit.spline$predDat

coef.Dat <- fit.spline$coefDat


out1 <- detectOutliers(corrDat = spatCorr,
                       predDat = pred.Dat,
                       coefDat = coef.Dat,
                       trait = "EffpsII",
                       genotypes =  c("G15","G121","G31","check2","check3","check4" ),# unique(spatCorr[["genotype"]]),
                       thrCor = 0.9,
                       thrPca = 1,
                       title = "Phenovator",
                       outFile = "testOutliers.pdf")

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
     whichED = c("repId:colId", "rowId", "fColRow","colfRow", "surface"),
     EDType = "ratio")
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
## plot Spat
plot(obj,
     timePoints = 41,
     plotType = "spatial",
     spaTrend = "raw") #"raw" #"percentage"
plot(obj,
     timePoints = 41,
     plotType = "spatial",
     spaTrend = "percentage") #"raw" #"percentage"


#----------------------------------------------------------------------------------------
#####-
##### Fitting ASReml models on the Phenovator data #######################################-
#####-

#####
##### ------ Fit a model few time points with geno = R ----------------------------------#
#####
modPhenoSpAs <- fitModels(TP = phenoTP,
                          trait = "pheno",
                          timePoints = seq(1,73,by=5),
                          engine = "asreml",
                          spatial = TRUE)

obj <- modPhenoSpAs
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 6)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 6)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
# genotypes = c("toto", "tutu", "tata"))
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit", yLim = c(0,1))
## plot Var
plot(obj,
     plotType = "variance")
# plot Spat
plot(obj,
     timePoints = 36,
     plotType = "spatial")


#####
##### ------ Fit a model few time points with geno = R and Covar ------------------------#
#####
modPhenoSpAsCov <- fitModels(TP = phenoTP,
                             trait = "pheno",
                             covariates = c("repId","Image_pos"),
                             timePoints = seq(1,73,by=10),
                             engine = "asreml",
                             spatial = TRUE)

obj <- modPhenoSpAsCov
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 11)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 11)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
# genotypes = c("toto", "tutu", "tata"))
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit", yLim = c(0,1))
## plot Var
plot(obj,
     plotType = "variance")
# plot Spat
plot(obj,
     timePoints = 36,
     plotType = "spatial")

#####
##### ------ Fit a model few time points with geno = R and check ------------------------#
#####
modPhenoSpAsCheck <- fitModels(TP = phenoTP,
                               trait = "pheno",
                               useCheck = TRUE,
                               timePoints = seq(1,73,by=10),
                               engine = "asreml",
                               spatial = TRUE)

obj <- modPhenoSpAsCheck
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 6)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 6)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("col", "ler", "ely", "evo1")) #levels(factor(PhenovatorDat1$Genotype))[1:28] )
# genotypes = c("toto", "tutu", "tata"))
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("col", "ler", "ely", "evo1") )
## plot Herit
plot(obj,
     plotType = "herit", yLim = c(0,1))
## plot Var
plot(obj,
     plotType = "variance")
# plot Spat
plot(obj,
     timePoints = 36,
     plotType = "spatial")





#####



#####-
##### Fitting ASReml models on the Phenoarch data #######################################-
#####-

data("PhenoarchDat1")
phenoTParch <- createTimePoints(dat = PhenoarchDat1,
                                experimentName = "ZA17",
                                genotype = "geno",
                                timePoint = "Date",
                                repId = "Rep",
                                plotId = "pos",
                                rowNum = "Position",
                                colNum = "Line")

modPhenoSpCov <- fitModels(TP = phenoTParch,
                           trait = "LA_Estimated",
                           # geno.decomp = c("Scenario","population"),
                           covariates = "Scenario",
                           timePoints = seq(1,35,by=10),
                           spatial = TRUE,
                           what = "fixed")

modPhenoSpAsGD <- fitModels(TP = phenoTParch,
                            trait = "LA_Estimated",
                            geno.decomp = c("Scenario","population"),
                            # covariates = "Image_pos",
                            timePoints = seq(1,35,by=10),
                            spatial = TRUE,
                            engine = "asreml",
                            what = "fixed")

obj <- modPhenoSpAsGD
obj <- modPhenoSpCov
# Extract the genotypic predictions:
genoPred <- getGenoPred(obj, timePoints = 11)
# Extract the corrected values:
spatCorr <- getCorrected(obj, timePoints = 11)
# Extract model components:
variance <- getVar(obj)
herit    <- getHerit(obj)
## plot Pred
plot(obj,
     plotType = "rawPred",
     genotypes = c("GenoA1","GenoA2","GenoB1","GenoB2") )
## plot Corr
plot(obj,
     plotType = "corrPred",
     genotypes = c("GenoA1","GenoA2","GenoB1","GenoB2")  )
## plot Herit
plot(obj,
     plotType = "herit", yLim = c(0,1))
## plot Var
plot(obj,
     plotType = "variance")
# plot Spat
plot(obj,
     timePoints = 11,
     plotType = "spatial")


#####-
##### Spline and outlier on Phenovator data #######################################-
#####-

modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = seq(from = 1, to = 73, by = 5))

# Extract the corrected values:
spatCorr <- getCorrected(modPhenoSp)

splineSp <- fitSpline(spatCorr,
                      "EffpsII_corr",
                      knots = 5)

