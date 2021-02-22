## Create data for vignette.

####### 1. Phenovator data set
#### 1.1 Raw data
# Read raw data.
PhenovatorDat1 <- read.csv("./data-raw/Phenovator_Data_Example1.csv",
                           stringsAsFactors = TRUE)
# Remove data with missing phenotype.
PhenovatorDat1 <- PhenovatorDat1[!is.na(PhenovatorDat1$EffpsII), ]
# Create an indicator for each plot (according to the row and column position).
PhenovatorDat1$pos <- paste0("c", PhenovatorDat1[["x"]],
                             "r", PhenovatorDat1[["y"]])
# Remove a plant that has very few measurements.
PhenovatorDat1 <- PhenovatorDat1[PhenovatorDat1$pos != "c1r54", ]
# Rename genotypes to G001 - G192
genoVator <- levels(PhenovatorDat1$Genotype)
genoVator[substr(genoVator, 1, 1) == "G"] <-
  paste0("G", formatC(as.numeric(substring(genoVator[substr(genoVator, 1, 1) == "G"], 2)),
                      digits = 2, flag = "0", format = "d"))
levels(PhenovatorDat1$Genotype) <- genoVator

# Export to package
usethis::use_data(PhenovatorDat1, overwrite = TRUE)

#### 1.2. Corrected data - outliers removed

## Create csv containing corrected data.
# Commented out since it runs a long time and the .csv is saved for later use.
# PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
#                                    c("c24r41", "c7r18", "c7r49"),]
# # Create TP object.
# phenoTP <- createTimePoints(dat = PhenovatorDat1, experimentName = "Phenovator",
#                             genotype = "Genotype", timePoint = "timepoints",
#                             repId = "Replicate", plotId = "pos", rowNum = "y",
#                             colNum = "x", addCheck = TRUE,
#                             checkGenotypes = c("check1", "check2",
#                                                "check3", "check4"))
# # Detect and remove outliers.
# resuVatorHTP <- detectSingleOut(TP = phenoTP, trait = "EffpsII",
#                                 confIntSize = 3, mylocfit = 0.1)
# phenoTPOut <- removeSingleOut(phenoTP, resuVatorHTP)
#
# # Fit models and get corrected values.
# modPhenoSpCheck <- fitModels(TP = phenoTPOut, trait = "EffpsII",
#                              extraFixedFactors = c("repId", "Image_pos"),
#                              useCheck = TRUE)
# spatCorrectedVator <- getCorrected(modPhenoSpCheck)
#
# # Write to .csv
# write.table(spatCorrectedVator, "./data-raw/PhenovatorDat1_corr_outPoint.csv",
#             sep = ",", row.names = FALSE)

# Read from .csv
spatCorrectedVator <- read.csv("./data-raw/PhenovatorDat1_corr_outPoint.csv",
                               stringsAsFactors = TRUE)
# Format the time in hour since first measurement
spatCorrectedVator$timePoint <- lubridate::as_datetime(spatCorrectedVator$timePoint)
timy <- data.frame(timePoint = unique(spatCorrectedVator$timePoint),
                   timePointP1 = c(unique(spatCorrectedVator$timePoint)[2:73],
                                   lubridate::ymd_hms("2018-06-18 16:37:00")),
                   timeNum = NA)
# diff between two time points in hours.
timeNum <- sapply(1:nrow(timy), function(x) {
  as.numeric(lubridate::ymd_hms(timy$timePointP1[x])-
               lubridate::ymd_hms(timy$timePoint)[x])
})
timy$timeNumDiff <- c(0,timeNum[1:(length(timeNum) - 1)])
# cum sum
timy$timeNum <- cumsum(timy$timeNumDiff)
# add to spatCorrVator
spatCorrectedVator$timeNumHour <- timy$timeNum[match(spatCorrectedVator$timePoint,
                                                     timy$timePoint)]
# Export to package
usethis::use_data(spatCorrectedVator, overwrite = TRUE)


####### 2. Phenoarch data set
#### 2.1 Raw data
# Read raw data.
PhenoarchDat1 <- read.csv("./data-raw/Phenoarch_ZA17_extraVariables.csv",
                          stringsAsFactors = TRUE)
# Rename genotypes.
genoArch <- levels(PhenoarchDat1$geno)
genoArchRen <- genoArch[nchar(genoArch) == 6]
genoArch[genoArch %in% genoArchRen] <-
  paste0(substr(genoArchRen, 1, 5), "0", substr(genoArchRen, 6, 6))
levels(PhenoarchDat1$geno) <- genoArch

# Export to package
usethis::use_data(PhenoarchDat1, overwrite = TRUE)

#### 2.2. Corrected data - outliers removed

## Create csv containing corrected data.
# Commented out since it runs a long time and the .csv is saved for later use.

# # Create TP object.
# phenoTParch <- createTimePoints(dat = PhenoarchDat1,
#                                 experimentName = "Phenoarch", genotype = "geno",
#                                 timePoint = "Time", plotId = "pos",
#                                 rowNum = "Row", colNum = "Col")
# # Detect and remove outliers.
# resuArchHTP <- detectSingleOut(TP = phenoTParch, trait = "LA_Estimated",
#                                confIntSize = 5, mylocfit = 0.5)
# # Fit models.
# phenoTParchOut <- removeSingleOut(phenoTParch, resuArchHTP)
#
# modPhenoSpGD <- fitModels(TP = phenoTParchOut, trait = "LA_Estimated",
#                           geno.decomp = c("Scenario", "population"))
#
# # Get and write corrected values and predictions.
# spatCorrectedArch <- getCorrected(modPhenoSpGD)
# write.table(spatCorrectedArch,
#             file = "./data-raw/PhenoArchDat1_corr_OutPoint_LA.csv",
#             sep = ",", row.names = FALSE)
# spatPredArch <- getGenoPred(modPhenoSpGD)$genoPred
# write.table(spatPredArch,
#             file = "./data-raw/PhenoArchDat1_pred_OutPoint_LA.csv",
#             sep = ",", row.names = FALSE)

# Read raw data.
spatCorrectedArch <- read.csv("./data-raw/PhenoArchDat1_corr_OutPoint_LA.csv",
                              stringsAsFactors = TRUE)
# Format the timepoint
spatCorrectedArch$timePoint <- lubridate::as_datetime(spatCorrectedArch$timePoint)
# Export to package
usethis::use_data(spatCorrectedArch, overwrite = TRUE)

#### 2.3. Genotypic prediction data - outliers removed
# Read raw data.
spatPredArch <- read.csv("./data-raw/PhenoArchDat1_pred_OutPoint_LA.csv",
                         stringsAsFactors = TRUE)
# Format the timepoint
spatPredArch$timePoint <- lubridate::as_datetime(spatPredArch$timePoint)
# Export to package
usethis::use_data(spatPredArch, overwrite = TRUE)


####### 3. RootUCL data set
#### 3.1 Raw data
# Read raw data.
RootDat1 <- read.csv("./data-raw/Data_tipclean_format.csv",
                     stringsAsFactors = TRUE)
# select one tank
RootDat1 <- RootDat1[RootDat1$Tank == "A", ]
# Export to package
usethis::use_data(RootDat1, overwrite = TRUE)

#### 3.2. not corrected data - outliers removed
# Read raw data.
noCorrectedRoot <- read.csv("./data-raw/RootDat1_nocorr.csv",
                       stringsAsFactors = TRUE)
noCorrectedRoot <- noCorrectedRoot[, c(11, 12, 3, 2, 4, 6, 5, 7, 8, 10)]
# Format the timepoint
noCorrectedRoot$timePoint <- lubridate::as_datetime(noCorrectedRoot$timePoint)
# Export to package
usethis::use_data(noCorrectedRoot, overwrite = TRUE)

## Create data for testing.

# Read raw data.
testDat <- read.csv("./data-raw/Phenovator_Data_Example1.csv",
                    stringsAsFactors = TRUE)
# Restrict data.
testDat <- with(testDat, testDat[x >= 10 & x < 15 & y <= 5 &
                                   as.numeric(timepoints) >= 5 &
                                   as.numeric(timepoints) <= 9, ])
# Rename trait column.
colnames(testDat)[colnames(testDat) == "EffpsII"] <- "t1"
# Create an indicator for each plot (according to the row and column position).
testDat$pos <- paste0("c", testDat[["x"]], "r", testDat[["y"]])
# Rename all checks to check1.
testDat$Genotype[testDat$Genotype %in%
                   c("check2", "check3", "check4", "G128")] <- "check1"
# Redefine replicates.
testDat$Replicate <- ifelse(testDat$pos %in%
                              c("c10r3", "c10r4", "c10r5", "c11r2", "c11r3",
                                "c11r4", "c11r5", "c12r2", "c12r3", "c12r4",
                                "c12r5", "c13r2"), 1, 2)
# Export to package
write.csv(testDat, file = "./inst/tinytest/testDat.csv", row.names = FALSE)
