### Test detectSingleOut.

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)

## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos")

## Check that general checks in detectSingleOut function correctly.
expect_error(detectSingleOut(1),
             "TP should be an object of class TP")
expect_error(detectSingleOut(testTP, trait = 1),
             "trait should be a character string of length 1")
expect_error(detectSingleOut(testTP, trait = c("t1", "t2")),
             "trait should be a character string of length 1")
expect_error(detectSingleOut(testTP, trait = "t2"),
             "TP should contain a column t2")
expect_error(detectSingleOut(testTP, trait = "t1", plotIds = "c1r1"),
             "All plotIds should be in TP")

## Check that check for minimal number of time points functions correctly.
expect_warning(detectSingleOut(testTP, trait = "t1"),
               "to fit a model for: c10r2")
expect_error(detectSingleOut(testTP, trait = "t1"),
             "for any of the plots")

## testDat only contains 5 timepoints so cannot be used for actual testing.
## Using the data in the package for that instead.
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            plotId = "pos")

## Select some plots.
plantSel <- c("c14r32", "c13r17")

singleOut <- detectSingleOut(phenoTP, trait = "EffpsII",
                             plotIds = plantSel)

## Check that general structure of the output is correct.
expect_inherits(singleOut, c("singleOut", "data.frame"))
expect_equal(dim(singleOut), c(132, 8))
expect_equal(colnames(singleOut),
             c("plotId", "timePoint", "EffpsII", "yPred", "sd_yPred",
               "lwr", "upr", "outlier"))

## Check that full output content is correct.
expect_equal_to_reference(singleOut, file = "singleOut", tolerance = 1e-5)

## Check that parameter confIntSize functions correctly.

# Setting a high value should lead to no outliers.
singleOut2 <- detectSingleOut(phenoTP, trait = "EffpsII",
                              plotIds = plantSel,
                              confIntSize = 15)
expect_equal(sum(singleOut2[["outlier"]]), 0)

# Setting confIntSize to 0 should cause (almost) all observations to be outliers.
singleOut3 <- detectSingleOut(phenoTP, trait = "EffpsII",
                              plotIds = plantSel,
                              confIntSize = 0)
expect_equal(sum(singleOut3[["outlier"]]), 132)

## Check that parameter nnLocFit functions correctly.

# A low value of nnLocFit should lead to a decrease of number of outliers.
singleOut4 <- detectSingleOut(phenoTP, trait = "EffpsII",
                              plotIds = plantSel,
                              nnLocfit = 0.1)
expect_equal(sum(singleOut4[["outlier"]]), 1)

## Check that parameter checkEdges functions correctly.

# Add an edge outlier.
phenoTP2 <- phenoTP
phenoTP2[[1]][phenoTP2[[1]][["plotId"]] == "c14r32", "EffpsII"] <- 1

# Not checking edges influences outliers detected elsewhere.
singleOut5 <- detectSingleOut(phenoTP2, trait = "EffpsII",
                              plotIds = "c14r32")
singleOut6 <- detectSingleOut(phenoTP2, trait = "EffpsII",
                              plotIds = "c14r32",
                              checkEdges = FALSE)
expect_equal(
  as.numeric(setdiff(singleOut5[singleOut5[["outlier"]] == 1, "timePoint"],
                     singleOut6[singleOut6[["outlier"]] == 1, "timePoint"])),
  1528448820)
expect_equal(
  as.numeric(setdiff(singleOut6[singleOut6[["outlier"]] == 1, "timePoint"],
                     singleOut5[singleOut5[["outlier"]] == 1, "timePoint"])),
  1527871020)

### Check plotting of detectSingleOut results.

## Check that general checks in plot function correctly.
expect_error(plot(singleOut, plotIds = "a"),
             "All plotIds should be in x")

## Check that general output structure is correct.
expect_silent(p <- plot(singleOut))
expect_inherits(p, "list")
expect_equal(length(p), 1)
expect_inherits(p[[1]], "ggplot")

## Check that parameter outOnly functions correctly.
# singleOut2 has no ouliers.

expect_error(plot(singleOut2),
             "No outliers present for selected plotIds")
expect_silent(plot(singleOut2, outOnly = FALSE))

### Check removal of outliers detected by detectSingleOut

## Check that general checks in plot function correctly.
expect_error(removeSingleOut(1),
             "TP should be an object of class TP")
expect_error(removeSingleOut(phenoTP, singleOut = 1),
             "singleOut should be a data.frame")
expect_error(removeSingleOut(phenoTP, singleOut = data.frame()),
             "singleOut should at least contain the columns plotId and timePoint")
expect_error(removeSingleOut(phenoTP[1], singleOut = singleOut),
             "All time points in singleOut should be in TP")
expect_error(removeSingleOut(phenoTP, singleOut = singleOut2),
             "There are no outlying points in TP")

## Check that outliers are removed.
phenoTPOut <- removeSingleOut(phenoTP, singleOut = singleOut)

expect_inherits(phenoTPOut, "TP")

expect_true(is.na(phenoTPOut[[11]][phenoTPOut[[11]][["plotId"]] == "c13r17",
                                   "EffpsII"]))

## Check that paramater trait functions correctly.
# Data contains only one trait, so set Sowing_Position to NA.
phenoTPOut2 <- removeSingleOut(phenoTP, singleOut = singleOut,
                               trait = "Sowing_Position")
expect_true(is.na(phenoTPOut2[[11]][phenoTPOut2[[11]][["plotId"]] == "c13r17",
                                    "Sowing_Position"]))

