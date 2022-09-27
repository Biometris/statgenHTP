### Test plot.psHDM

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
## Move all check1s to Basin 1 for a consistent population structure.
corr[corr[["genotype"]] == "check1", "Basin"] <- 1

## Fit a HDM spline on the corrected values.
splineRes <- fitSplineHDM(inDat = corr, trait = "t1_corr", pop = "Basin",
                          trace = FALSE)

## Make most restrictive and broadest predictions.
pred <- predict(splineRes,
                pred = list(pop = TRUE, geno = FALSE, plot = FALSE),
                se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                trace = FALSE)
pred1 <- predict(splineRes,
                 pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                 se = list(pop = TRUE, geno = TRUE, plot = TRUE),
                 trace = FALSE)


## Create tempfile for plotting output.
tmpFile <- tempfile(fileext = ".pdf")

## General checks.
expect_error(plot(splineRes, xlab = 1:2),
             "xlab should have length 1")
expect_error(plot(splineRes, ylab = 1:2),
             "ylab should have length 1")
expect_error(plot(splineRes, title = 1:2),
             "title should have length 1")

expect_error(plot(pred, plotType = "popGenoTra"),
             paste("Genotype-specific growth curves can only be plotted if",
                   "predictions were made at genotype and population level"))
expect_error(plot(pred, plotType = "popGenoDeriv"),
             paste("First order derivatives of genotype-specific growth curves",
                   "can only be plotted if predictions were made at genotype level"))
expect_error(plot(pred, plotType = "genoDev"),
             paste("Genotype-specific deviations can only be plotted if",
                   "predictions were made at genotype level"))
expect_error(plot(pred, plotType = "genoPlotTra"),
             paste("Plot and Genotype-specific growth curves can only be",
                   "plotted if predictions were made at genotype and plot level"))


## Check output.
p <- plot(pred1, plotType = "popTra", outFile = tmpFile)
expect_inherits(p, "ggplot")

p2 <- plot(pred1, plotType = "popGenoTra", outFile = tmpFile)
expect_inherits(p2, "ggplot")

p3 <- plot(pred1, plotType = "popGenoDeriv", outFile = tmpFile)
expect_inherits(p3, "ggplot")

p4 <- plot(pred1, plotType = "genoDev", outFile = tmpFile)
expect_inherits(p4, "ggplot")

p5 <- plot(pred1, plotType = "genoPlotTra", outFile = tmpFile)
expect_inherits(p5, "list")
expect_inherits(p5[[1]], "ggplot")


## Check that option genotypes functions correctly.
expect_error(plot(splineRes, plotType = "genoPlotTra", genotypes = "G0"),
             "All genotypes should be in fitted model")
p6 <- plot(splineRes, plotType = "genoPlotTra", genotypes = "check1",
           outFile = tmpFile)
expect_equal(nrow(p6[[1]]$data), 20)

## Check that option genotypeNames functions correctly.
p7 <- plot(splineRes, plotType = "genoPlotTra", genotypes = "check1",
           genotypeNames = "a", outFile = tmpFile)
expect_equal(levels(p7[[1]]$data$genotype), "a")

## Check that option genotypeOrder functions correctly.
p8 <- plot(splineRes, plotType = "genoPlotTra", genotypes = c("check1", "G12"),
           genotypeOrder = 2:1, outFile = tmpFile)
expect_equal(levels(p8[[1]]$data$genotype), c("check1", "G12"))
