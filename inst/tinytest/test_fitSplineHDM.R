### Test fitSpline

Sys.setlocale("LC_COLLATE", "C")

## Read test data from .csv
testDat <- read.csv("testDat.csv", stringsAsFactors = FALSE)

## Create TP object.
testTP <- createTimePoints(dat = testDat, experimentName = "testExp",
                           genotype = "Genotype", timePoint = "timepoints",
                           repId = "Replicate", plotId = "pos",
                           rowNum = "y", colNum = "x")

## Fit model.
testFitMod <- fitModels(testTP, trait = "t1", geno.decomp = "repId",
                        extraFixedFactors = "Basin", quiet = TRUE)

## Get corrected values.
corr <- getCorrected(testFitMod)

## Check that general checks in fitSpline function correctly.
expect_error(fitSplineHDM(1, trait = "t1_corr"),
             "inDat should be a data.frame")
expect_error(fitSplineHDM(corr, trait = 1),
             "trait should be a character string of length 1")
expect_error(fitSplineHDM(corr, trait = "t1_corr"),
             "inDat should at least contain the following columns")
expect_error(fitSplineHDM(corr, trait = "t1_corr", pop = "Basin"),
             "The following genotypes are in multiple populations")

## Move all check1s to Basin 1 for a consistent population structure.
corr[corr[["genotype"]] == "check1", "Basin"] <- 1

## Move 1 plot to the wrong genotype.
corrErr <- corr
corrErr[which(corrErr[["plotId"]] == "c10r3")[1], "genotype"] <- "G3"

expect_error(fitSplineHDM(corrErr, trait = "t1_corr", pop = "Basin"),
             "The following plots are specified for multiple genotypes")


## Fit a HDM spline on the corrected values.
expect_warning(splineRes <- fitSplineHDM(inDat = corr, trait = "t1_corr",
                                         pop = "Basin", trace = FALSE),
               "for less than the minimum number of time points, which is 3")

## Check that general structure of the output is correct.
expect_inherits(splineRes, "psHDM")
expect_equal(length(splineRes), 22)
expect_equal(names(splineRes),
             c("y", "time", "popLevs", "genoLevs", "plotLevs", "nPlotPop",
               "nGenoPop", "nPlotGeno", "MM", "ed", "vc", "phi", "coeff",
               "deviance", "convergence", "dim", "family", "Vp", "smooth",
               "popLevel", "genoLevel", "plotLevel"))

## Check structure of components.
expect_inherits(splineRes[["y"]], "list")
expect_equal(length(splineRes[["y"]]), 22)

expect_inherits(splineRes[["time"]], "data.frame")
expect_inherits(splineRes[["genoLevs"]], "factor")
expect_inherits(splineRes[["popLevs"]], "factor")
expect_inherits(splineRes[["plotLevs"]], "factor")
expect_inherits(splineRes[["nPlotPop"]], "integer")
expect_inherits(splineRes[["nGenoPop"]], "integer")
expect_inherits(splineRes[["nPlotGeno"]], "integer")

expect_inherits(splineRes[["MM"]], "list")
expect_equal(length(splineRes[["MM"]]), 3)
expect_equal(names(splineRes[["MM"]]),
             c("MMPop", "MMGeno", "MMPlot"))

expect_inherits(splineRes[["ed"]], "numeric")
expect_equal(length(splineRes[["ed"]]), 8)

expect_inherits(splineRes[["vc"]], "numeric")
expect_equal(length(splineRes[["vc"]]), 8)

expect_inherits(splineRes[["phi"]], "numeric")

expect_inherits(splineRes[["coeff"]], "numeric")
expect_equal(length(splineRes[["coeff"]]), 637)

expect_inherits(splineRes[["deviance"]], "numeric")
expect_inherits(splineRes[["convergence"]], "logical")

expect_inherits(splineRes[["dim"]], "integer")
expect_equal(length(splineRes[["dim"]]), 6)

expect_inherits(splineRes[["family"]], "family")

expect_inherits(splineRes[["Vp"]], "matrix")
expect_equal(dim(splineRes[["Vp"]]), c(637, 637))

expect_inherits(splineRes[["smooth"]], "list")
expect_equal(length(splineRes[["smooth"]]), 3)
expect_equal(names(splineRes[["smooth"]]),
             c("smoothPop", "smoothGeno", "smoothPlot"))

expect_inherits(splineRes[["popLevel"]], "data.frame")
expect_inherits(splineRes[["genoLevel"]], "data.frame")
expect_inherits(splineRes[["plotLevel"]], "data.frame")


## get predictions at all levels.
popLevel <- splineRes[["popLevel"]]
genoLevel <- splineRes[["genoLevel"]]
plotLevel <- splineRes[["plotLevel"]]

## Check consistency with other outputs.
expect_equal(levels(popLevel[["pop"]]), levels(splineRes[["popLevs"]]))

expect_equal(levels(genoLevel[["pop"]]), levels(splineRes[["popLevs"]]))
expect_equal(levels(genoLevel[["genotype"]]), levels(splineRes[["genoLevs"]]))

expect_equal(levels(plotLevel[["pop"]]), levels(splineRes[["popLevs"]]))
expect_equal(levels(plotLevel[["genotype"]]), levels(splineRes[["genoLevs"]]))
expect_equal(levels(plotLevel[["plotId"]]), levels(splineRes[["plotLevs"]]))

## Check that full prediction results are correct.
expect_equal_to_reference(popLevel, file = "popLevel", tolerance = 1e-6)
expect_equal_to_reference(genoLevel, file = "genoLevel", tolerance = 1e-6)
expect_equal_to_reference(plotLevel, file = "plotLevel", tolerance = 1e-6)


## Check that option timeNumber functions correctly.
expect_error(fitSplineHDM(corr, trait = "t1_corr", pop = "Basin",
                          useTimeNumber = TRUE),
             "timeNumber should be a character string of length 1")
expect_error(fitSplineHDM(corr, trait = "t1_corr", pop = "Basin",
                          useTimeNumber = TRUE, timeNumber = "genotype"),
             "timeNumber should be a numerical column")

expect_warning(fitSplineHDM(inDat = corr, trait = "t1_corr",
                            pop = "Basin", useTimeNumber = TRUE,
                            timeNumber = "timeNumber", trace = FALSE),
               "for less than the minimum number of time points, which is 3")


## Check that option genotypes functions correctly.
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          genotypes = 1),
             "genotypes should be a character vector of genotypes in inDat")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          genotypes = "G0"),
             "genotypes should be a character vector of genotypes in inDat")

splineRes2 <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                           genotypes = "G12", trace = FALSE)

expect_true(all(splineRes2[["genoLevel"]][["genotype"]] == "G12"))
expect_true(all(splineRes2[["plotLevel"]][["genotype"]] == "G12"))



## Check that option plotIds functions correctly.
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          plotIds = 1),
             "plotIds should be a character vector of plotIds in inDat")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          plotIds = "p1"),
             "plotIds should be a character vector of plotIds in inDat")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          genotypes = "G12", plotIds = "c13r2"),
             "At least one valid combination of genotype and plotId should be selected")

splineRes3 <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                           plotIds = "c13r2", trace = FALSE)

expect_true(all(splineRes3[["plotLevel"]][["plotId"]] == "c13r2"))


## Check that option minNoTP functions correctly.
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          minNoTP = "a"),
             "minNoTP should be a numerical value")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          minNoTP = 50),
             "minNoTP should be a number bewtween 0 and 5")
expect_silent(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                           genotypes = "G12", minNoTP = 5))

# Add some extra NA to corr.
corr2 <- corr
corr2[1:10, "t1_corr"] <- NA

expect_warning(splineRes4 <- fitSplineHDM(inDat = corr2, trait = "t1_corr",
                                          pop = "Basin", minNoTP = 5,
                                          trace = FALSE),
               "More than 5 plots have observations for less than the minimum")
expect_equal(attr(splineRes4, which = "plotLimObs"),
             c("c10r1", "c10r3", "c10r4", "c10r5", "c11r2", "c11r3", "c11r4",
               "c11r5", "c12r3", "c12r4", "c12r5"))


## Check that option maxit functions correctly.
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          maxit = "a"),
             "maxit should be a positive numerical value")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          maxit = -1),
             "maxit should be a positive numerical value")

splineRes5 <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                           maxit = 1, trace = FALSE)
expect_false(splineRes5[["convergence"]])


## Check that option thr functions correctly.
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          thr = "a"),
             "thr should be a positive numerical value")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          thr = -1),
             "thr should be a positive numerical value")

splineRes6 <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                           thr = 1, trace = FALSE)
expect_true(splineRes6[["deviance"]] > splineRes[["deviance"]])



## Check that difVar functions correctly.

expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr",
                          pop = "Basin", difVar = "a"),
             "difVar should be a named list of length 2")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr",
                          pop = "Basin", difVar = list(geno = TRUE)),
             "difVar should be a named list of length 2")
expect_error(fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          difVar = list(geno = TRUE, pop = FALSE)),
             "difVar should be a named list of length 2")

splineRes7a <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                            difVar = list(geno = TRUE, plot = FALSE),
                            trace = FALSE)
splineRes7b <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                            difVar = list(geno = FALSE, plot = TRUE),
                            trace = FALSE)
splineRes7c <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                            difVar = list(geno = TRUE, plot = TRUE),
                            trace = FALSE)

expect_equal(length(splineRes7a[["vc"]]), 11)
expect_equal(length(splineRes7b[["vc"]]), 71)
expect_equal(length(splineRes7c[["vc"]]), 74)

