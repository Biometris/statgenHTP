## Create data for vignette.

####### 1. Phenovator data set
#### 1.1 Raw data
# Read raw data.
PhenovatorDat1 <- read.csv(system.file("extdata", "Phenovator_Data_Example1.csv",
                                       package = "statgenHTP"))
# Remove data with missing phenotype.
PhenovatorDat1 <- PhenovatorDat1[!is.na(PhenovatorDat1$EffpsII), ]
# Create an indicator for each plot (according to the row and column position).
PhenovatorDat1$pos <- paste0("c", PhenovatorDat1[["x"]],
                             "r", PhenovatorDat1[["y"]])
# Remove a plant that has very few measurements.
PhenovatorDat1 <- PhenovatorDat1[PhenovatorDat1$pos != "c1r54",]
# Export to package
usethis::use_data(PhenovatorDat1, overwrite = TRUE)

#### 1.2. Corrected data - outliers removed
# Read raw data.
spatCorrVator <- read.csv(system.file("extdata", "PhenovatorDat1_corr_outPoint.csv",
                                       package = "statgenHTP"))
# Format the time in hour since first measurement
spatCorrVator$timePoint <- lubridate::as_datetime(spatCorrVator$timePoint)
timy <- data.frame(timePoint = unique(spatCorrVator$timePoint),
                   timePointP1 = c(unique(spatCorrVator$timePoint)[2:73],
                                   lubridate::ymd_hms("2018-06-18 16:37:00")),
                   timeNum = NA)
# diff between two time point in hour
timeNum <- sapply( 1:nrow(timy), function(x) {
  as.numeric(lubridate::ymd_hms(timy$timePointP1[x])-
               lubridate::ymd_hms(timy$timePoint)[x])
})
timy$timeNumDiff <- c(0,timeNum[1:(length(timeNum)-1)])
# cum sum
timy$timeNum <- cumsum(timy$timeNumDiff)
# add to spatCorrVator
spatCorrVator$timeNumHour <- timy$timeNum[match(spatCorrVator$timePoint,timy$timePoint)]
# Export to package
usethis::use_data(spatCorrVator, overwrite = TRUE)


####### 2. Phenoarch data set
#### 2.1 Raw data
# Read raw data.
# PhenoarchDat1 <- read.csv(system.file("extdata", "Phenoarch_Data_ZA17.csv",
#                                       package = "statgenHTP"))
PhenoarchDat1 <- read.csv(system.file("extdata", "Phenoarch_ZA17_extraVariables.csv",
                                      package = "statgenHTP"))
# Create an indicator for each plot (according to the row and column position)
# PhenoarchDat1$pos <- paste0("c", PhenoarchDat1$Line, "r", PhenoarchDat1$Position)
# Export to package
usethis::use_data(PhenoarchDat1, overwrite = TRUE)


#### 2.2. Corrected data - outliers removed
# Read raw data.
spatCorrArch <- read.csv(system.file("extdata", "PhenoarchDat1_corr_outPoint.csv",
                                      package = "statgenHTP"))
# Export to package
usethis::use_data(spatCorrArch, overwrite = TRUE)


####### 3. RootUCL data set
#### 3.1 Raw data
# Read raw data.
RootDat1 <- read.csv(system.file("extdata", "Data_tipclean_format.csv",
                                      package = "statgenHTP"))
# select one tank
RootDat1 <- RootDat1[RootDat1$Tank == "A",]
# Export to package
usethis::use_data(RootDat1, overwrite = TRUE)

#### 3.2. not corrected data - outliers removed
# Read raw data.
noCorrRoot <- read.csv(system.file("extdata", "RootDat1_nocorr.csv",
                                     package = "statgenHTP"))
noCorrRoot <- noCorrRoot[,c(11,12,3,2,4,6,5,7,8,10)]
# Export to package
usethis::use_data(noCorrRoot, overwrite = TRUE)



## Create data for testing.

# Read raw data.
testDat <- read.csv(system.file("extdata",
                                "Phenovator_Data_Example1.csv",
                                package = "statgenHTP"))
# Restrict data.
testDat <- with(testDat, testDat[x >= 10 & x < 15 & y <= 5 &
                                   as.numeric(timepoints) >= 5 &
                                   as.numeric(timepoints) <= 9, ])
# Rename trait column.
colnames(testDat)[colnames(testDat) == "EffpsII"] <- "t1"
# Create an indicator for each plot (according to the row and column position).
testDat$pos <- paste0("c", testDat[["x"]], "r", testDat[["y"]])
# Export to package
write.csv(testDat, file = "./inst/tinytest/testDat.csv", row.names = FALSE)
