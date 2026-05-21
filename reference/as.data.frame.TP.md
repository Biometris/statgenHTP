# Coerce TP object to data.frame

Function for converting an object of class TP to a data.frame.

## Usage

``` r
# S3 method for class 'TP'
as.data.frame(x, ...)
```

## Arguments

- x:

  An object of class TP.

- ...:

  Ignored.

## Value

A data.frame containing the data.frames for all time points in the TP
object bound together.

## See also

Other functions for data preparation:
[`createTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/createTimePoints.md),
[`getTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/getTimePoints.md),
[`plot.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.TP.md),
[`removeTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/removeTimePoints.md),
[`summary.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.TP.md)

## Examples

``` r
## Create a TP object containing the data from the Phenovator.
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
#> Warning: The following plotIds have observations for less than 50% of the time points:
#> c24r41, c7r18, c7r49
## Convert phenoTP to data.frame.
phenoDat <- as.data.frame(phenoTP)
```
