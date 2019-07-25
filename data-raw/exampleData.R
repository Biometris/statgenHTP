## Create data for vignette.

# Read raw data.
PhenovatorDat1 <- read.csv(system.file("extdata", "Original_PAM_reshape.csv",
                                       package = "statgenHTP"))
# Remove data with missing phenotype.
PhenovatorDat1 <- PhenovatorDat1[!is.na(PhenovatorDat1$pheno), ]
# Create an indicator for each plot (according to the row and column position).
PhenovatorDat1$pos <- paste0("c", PhenovatorDat1[["x"]],
                             "r", PhenovatorDat1[["y"]])
# Remove a plant that has very few measurements.
PhenovatorDat1 <- PhenovatorDat1[PhenovatorDat1$pos != "c1r54",]

# Export to package
usethis::use_data(PhenovatorDat1, overwrite = TRUE)
