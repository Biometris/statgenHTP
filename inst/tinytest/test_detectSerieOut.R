### Test detectSerieOut.

Sys.setlocale("LC_COLLATE", "C")

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)
# Set Replicate to 1 for all observations of check1 so geno.decomp can
# actually be used for testing.
testDat[testDat[["Genotype"]] == "check1", "Replicate"] <- 1

## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           plotId = "pos", rowNum = "y", colNum = "x")

## Fit models.
testFitMod <- fitModels(testTP, trait = "t1", quiet = TRUE)
testFitModGD <- fitModels(testTP, trait = "t1", geno.decomp = "Replicate",
                          quiet = TRUE)

## Get corrected values.
corr <- getCorrected(testFitMod)
corrGD <- getCorrected(testFitModGD)

## Fit a splines on the corrected values.
expect_warning(splineRes <- fitSpline(inDat = corr, trait = "t1_corr"))
expect_warning(splineResGD <- fitSpline(inDat = corrGD, trait = "t1_corr"))

## get coefDat and predDat.
coefDat <- splineRes[["coefDat"]]
predDat <- splineRes[["predDat"]]

coefDatGD <- splineResGD[["coefDat"]]
predDatGD <- splineResGD[["predDat"]]

## Check that general checks in detectSingleOut function correctly.
expect_error(detectSerieOut(trait = 1),
             "trait should be a character string of length 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = 1),
             "corrDat should be a data.frame")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr, predDat = 1),
             "predDat should be a data.frame")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = 1),
             "coefDat should be a data.frame")
expect_error(detectSerieOut(trait = "t2_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat),
             "corrDat should at least contain the following columns")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = coefDat, coefDat = coefDat),
             "predDat should at least contain the following columns")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = predDat),
             "coefDat should at least contain the following columns")

## Check that checks for genotypes function correctly
expect_warning(detectSerieOut(trait = "t1_corr", corrDat = corr,
                              predDat = predDat, coefDat = coefDat),
               "have less than 3 plotIds and are skipped in the outlier")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = 1),
             "genotypes should be a character vector")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "a"),
             "all genotypes should be in predDat")
expect_error(detectSerieOut(trait = "t1_corr",
                            corrDat = corr[corr[["genotype"]] != "G12", ],
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "G12"),
             "all genotypes should be in corrDat")
expect_error(detectSerieOut(trait = "t1_corr",  corrDat = corr,
                            predDat = predDat,
                            coefDat = coefDat[coefDat[["genotype"]] != "G12", ],
                            genotypes = "G12"),
             "all genotypes should be in coefDat")
expect_error(detectSerieOut(trait = "t1_corr",  corrDat = corr,
                            predDat = predDat,  coefDat = coefDat,
                            genotypes = "G12"),
             "No genotypes left for performing outlier detection")
expect_silent(serieOut1 <- detectSerieOut(trait = "t1_corr", corrDat = corr,
                                          predDat = predDat, coefDat = coefDat,
                                          genotypes = "check1"))

## Check that general structure of the output is correct.
expect_inherits(serieOut1, c("serieOut", "data.frame"))
expect_equal(dim(serieOut1), c(7, 4))
expect_equal(colnames(serieOut1),
             c("plotId", "genotype", "reason", "value"))

## Check that full output content is correct.
expect_equal_to_reference(serieOut1, file = "serieOut", tolerance = 1e-5)

## Check that parameter thrCor functions correctly.
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrCor = "a"),
             "thrCor should be a numerical vector with values between -1 and 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrCor = 2),
             "thrCor should be a numerical vector with values between -1 and 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrCor = 0:1),
             "thrCor should be a vector of length 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corrGD,
                            predDat = predDatGD, coefDat = coefDatGD,
                            genotypes = "check1", geno.decomp = "geno.decomp",
                            thrCor = 0:1),
             "thrCor should be a named vector, with names matching the levels")

## Check that parameter thrPca functions correctly.
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrPca = "a"),
             "thrPca should be a numerical vector with positive values")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrPca = -1),
             "thrPca should be a numerical vector with positive values")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrPca = 0:1),
             "thrPca should be a vector of length 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corrGD,
                            predDat = predDatGD, coefDat = coefDatGD,
                            genotypes = "check1", geno.decomp = "geno.decomp",
                            thrPca = 0:1),
             "thrPca should be a named vector, with names matching the levels")


## Check that parameter thrPca functions correctly.
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrSlope = "a"),
             "thrSlope should be a numerical vector with values between 0 and 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrSlope = -1),
             "thrSlope should be a numerical vector with values between 0 and 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrSlope = 0:1),
             "thrSlope should be a vector of length 1")
expect_error(detectSerieOut(trait = "t1_corr", corrDat = corrGD,
                            predDat = predDatGD, coefDat = coefDatGD,
                            genotypes = "check1", geno.decomp = "geno.decomp",
                            thrSlope = 0:1),
             "thrSlope should be a named vector, with names matching the levels")

# Set all three to extreme values should result in either 0 or many outliers.
serieOut2 <- detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrCor = 1, thrPca = 0,
                            thrSlope = 1)
expect_equal(dim(serieOut2), c(9, 4))

serieOut3 <- detectSerieOut(trait = "t1_corr", corrDat = corr,
                            predDat = predDat, coefDat = coefDat,
                            genotypes = "check1", thrCor = -1, thrPca = 180,
                            thrSlope = 0)
expect_equivalent(serieOut3, data.frame())

## Check that detecting outliers functions correctly for models with geno.decomp
expect_silent(serieOutGD <-
                detectSerieOut(trait = "t1_corr", corrDat = corrGD,
                               predDat = predDatGD, coefDat = coefDatGD,
                               genotypes = "check1", geno.decomp = "geno.decomp"))
expect_equal(dim(serieOutGD), c(7, 5))

## Check detectSerieOut functions correctly when plotIds are numeric-like.

corrNw <- corr
plotIds <- unique(corrNw$plotId)
corrNw[["plotId"]] <- as.factor(match(x = corrNw[["plotId"]], table = plotIds))
predDatNw <- predDat
predDatNw[["plotId"]] <- as.character(match(x = predDatNw[["plotId"]], table = plotIds))
coefDatNw <- coefDat
coefDatNw[["plotId"]] <- as.character(match(x = coefDatNw[["plotId"]], table = plotIds))

expect_silent(serieOut4 <- detectSerieOut(trait = "t1_corr", corrDat = corrNw,
                                          predDat = predDatNw,
                                          coefDat = coefDatNw,
                                          genotypes = "check1"))

## Variables in cormats and slopemats should be converted to factors for plotting.
expect_inherits(attr(serieOut4, which = "cormats")[[1]][["Var1"]],
                "factor")
expect_inherits(attr(serieOut4, which = "cormats")[[1]][["Var2"]],
                "factor")

expect_inherits(attr(serieOut4, which = "slopemats")[[1]][["Var1"]],
                "factor")
expect_inherits(attr(serieOut4, which = "slopemats")[[1]][["Var2"]],
                "factor")

### Check plotting of detectSerieOut results.

## Check that general checks in plot function correctly.
expect_error(plot(serieOut1, genotypes = "a"),
             "a character vector of genotypes used for outlier detection")

## Check that general output structure is correct.
expect_silent(p <- plot(serieOut1))
expect_inherits(p, "list")
expect_equal(length(p), 1)
expect_inherits(p[[1]], "gtable")

## Check that option useTimeNumber functions correctly.
expect_error(plot(serieOut1, useTimeNumber = TRUE, timeNumber = 1),
             "timeNumber should be a character string of length 1")
expect_silent(plot(serieOut1, useTimeNumber = TRUE, timeNumber = "timeNumber"))

## Check that option title functions correctly.
expect_silent(plot(serieOut1, genotypes = "check1", title = "bla"))

## Check that plotting functions correctly with geno.decomp.
expect_error(plot(serieOutGD, genotypes = "check1", geno.decomp = "2"),
             "All selected geno.decomp levels should be in the data")
expect_silent(plot(serieOutGD, genotypes = "check1", geno.decomp = "1"))

## Check that option reason functions correctly.
expect_error(plot(serieOut1, reason = "tst"),
             'one of "mean corr", "angle", "slope"')
expect_silent(plot(serieOut1, reason = "slope"))
expect_silent(plot(serieOut1, reason = c("slope", "angle")))

### Check removal of outliers detected by detectSerieOut

## Check that general checks in plot function correctly.
expect_error(removeSerieOut(dat = 1),
             "dat should be a data.frame")
expect_error(removeSerieOut(fitSpline = 1),
             "fitSpline should be an object of class HTPSpline")
expect_error(removeSerieOut(dat = corr, fitSpline = splineRes),
             "Specify exactly one of dat and fitSpline as inputs")
expect_error(removeSerieOut(dat = corr[, colnames(corr) != "plotId"]),
             "dat should at least contain the column plotId")
expect_error(removeSerieOut(dat = corr, serieOut = 1),
             "serieOut should be a data.frame")
expect_error(removeSerieOut(dat = corr,
                            serieOut = serieOut1[, colnames(serieOut1) != "plotId"]),
             "serieOut should at least contain the column plotId")

## Check that outliers are removed from data.frame.
corrOut1 <- removeSerieOut(dat = corr, serieOut = serieOut1)
expect_true(all(is.na(corrOut1[corrOut1[["plotId"]] == "c12r1", "t1_corr"])))

corrOut2 <- removeSerieOut(dat = corr, serieOut = serieOut1,
                           traits = c("t1", "t1_corr"))
expect_true(all(is.na(corrOut2[corrOut2[["plotId"]] == "c12r1",
                               c("t1", "t1_corr")])))

## Check that outliers are removed from fitted spline.
corrOut3 <- removeSerieOut(fitSpline = splineRes, serieOut = serieOut1)
coefOut3 <- corrOut3$coefDat
predOut3 <- corrOut3$predDat
modOut3 <- attr(corrOut3, "modDat")

expect_true(all(is.na(coefOut3[coefOut3[["plotId"]] == "c12r1",
                               "obj.coefficients"])))
expect_true(all(is.na(predOut3[predOut3[["plotId"]] == "c12r1",
                               c("pred.value", "deriv", "deriv2")])))
expect_true(all(is.na(modOut3[modOut3[["plotId"]] == "c12r1", "t1_corr"])))

## Check the option reason works correctly.
expect_error(removeSerieOut(dat = corr, serieOut = serieOut1, reason = "tst"),
             'one of "mean corr", "angle", "slope"')
expect_error(removeSerieOut(dat = corr,
                            serieOut = serieOut1[, colnames(serieOut1) != "reason"],
                            reason = "slope"),
             "serieOut should contain a column reason")
corrOut4 <- removeSerieOut(dat = corr, serieOut = serieOut1, reason = "angle")
expect_true(all(is.na(corrOut4[corrOut4[["plotId"]] == "c12r1", "t1_corr"])))

corrOut5 <- removeSerieOut(dat = corr, serieOut = serieOut1, reason = "slope")
expect_true(all(is.na(corrOut5[corrOut5[["plotId"]] == "c12r1", "t1_corr"])))

