# Remove time points from an object of class TP

Function for removing selected time points from an object of class TP.

## Usage

``` r
removeTimePoints(TP, timePoints)
```

## Arguments

- TP:

  An object of class TP.

- timePoints:

  A character or numeric vector indicating the time points to be
  removed. When using a character string to reference a time point, the
  value has to be an exact match to one of the existing timePoints. When
  using a number it will be matched by its number ("timeNumber") in the
  timePoints attribute of the TP object.

## Value

An object of class TP, the input with the selected time points removed.

## See also

Other functions for data preparation:
[`as.data.frame.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/as.data.frame.TP.md),
[`createTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/createTimePoints.md),
[`getTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/getTimePoints.md),
[`plot.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.TP.md),
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
                                               "check3","check4"))
#> Warning: The following plotIds have observations for less than 50% of the time points:
#> c24r41, c7r18, c7r49
## Remove the first and last time point from the TP object.
phenoTPNew <- removeTimePoints(phenoTP,
                               timePoints = c(1, 73))

## Compare by looking at summaries.
summary(phenoTP)
#> phenoTP contains data for experiment Phenovator.
#> 
#> It contains 73 time points.
#> First time point: 2018-05-31 16:37:00 
#> Last time point: 2018-06-18 16:37:00 
#> 
#> The following genotypes are defined as check genotypes: check1, check2, check3, check4.
summary(phenoTPNew)
#> phenoTPNew contains data for experiment Phenovator.
#> 
#> It contains 71 time points.
#> First time point: 2018-06-01 09:07:00 
#> Last time point: 2018-06-18 14:37:00 
#> 
#> No check genotypes are defined.
```
