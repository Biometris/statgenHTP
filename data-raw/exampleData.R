## Create data for vignette.

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


# Read raw data.
PhenoarchDat1 <- read.csv(system.file("extdata", "Phenoarch_Data_ZA17.csv",
                                      package = "statgenHTP"))
# Create an indicator for each plot (according to the row and column position)
PhenoarchDat1$pos <- paste0("c", PhenoarchDat1$Line, "r", PhenoarchDat1$Position)

# Export to package
usethis::use_data(PhenoarchDat1, overwrite = TRUE)


## Create data for testing.

# Read raw data.
testDat <- read.csv(system.file("extdata",
                                "Phenovator_Data_Example1.csv",
                                package = "statgenHTP"))
# Restrict data.
testDat <- with(testDat, testDat[x <= 5 & y <= 5 &
                                   as.numeric(timepoints) >= 5 &
                                   as.numeric(timepoints) <= 9, ])
# Create an indicator for each plot (according to the row and column position).
testDat$pos <- paste0("c", testDat[["x"]], "r", testDat[["y"]])
# Export to package
usethis::use_data(testDat, internal = TRUE, overwrite = TRUE)
