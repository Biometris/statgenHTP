# Plot outliers for single observations

Plot the fitted local regression, confidence intervals and detected
outliers for each plotId.

## Usage

``` r
# S3 method for class 'singleOut'
plot(x, ..., plotIds = NULL, outOnly = TRUE, output = TRUE)
```

## Arguments

- x:

  An object of class singleOut.

- ...:

  Ignored.

- plotIds:

  A character vector of plotIds for which the outliers should be
  detected. If `NULL`, all plotIds in `TP` are used.

- outOnly:

  Should only plots containing outliers be plotted?

- output:

  Should the plot be output to the current device? If `FALSE` only a
  (list of) ggplot object(s) is invisibly returned. Ignored if `outFile`
  is specified.

## Value

A list of ggplot objects is invisibly returned.

## See also

Other functions for detecting outliers for single observations:
[`detectSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOut.md),
[`detectSingleOutMaize()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOutMaize.md),
[`removeSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/removeSingleOut.md)

## Examples

``` r
## Create a TP object containing the data from the Phenovator.
PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
                                 c("c24r41", "c7r18", "c7r49"), ]
phenoTP <- createTimePoints(dat = PhenovatorDat1,
                            experimentName = "Phenovator",
                            genotype = "Genotype",
                            timePoint = "timepoints",
                            repId = "Replicate",
                            plotId = "pos",
                            rowNum = "y", colNum = "x",
                            addCheck = TRUE,
                            checkGenotypes = c("check1", "check2",
                                               "check3", "check4"))

## Select a subset of plants, for example here 9 plants.
plantSel <- phenoTP[[1]]$plotId[1:9]
# Then run on the subset.
resuVatorHTP <- detectSingleOut(TP = phenoTP,
                               trait = "EffpsII",
                               plotIds = plantSel,
                               confIntSize = 3,
                               nnLocfit = 0.1)

## Visualize the prediction by choosing a single plant...
plot(resuVatorHTP, plotIds = "c21r24", outOnly = FALSE)

## ...or a subset of plants.
plot(resuVatorHTP, plotIds = plantSel, outOnly = FALSE)

```
