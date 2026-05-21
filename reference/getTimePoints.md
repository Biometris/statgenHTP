# Extract time points

Function for extracting a data.frame with timeNumbers and timePoints
from an object of class TP or fitMod.

## Usage

``` r
getTimePoints(x)
```

## Arguments

- x:

  An object of class TP or fitMod

## Value

A data.frame with columns timeNumber and timePoint listing the time
points in x

## See also

Other functions for data preparation:
[`as.data.frame.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/as.data.frame.TP.md),
[`createTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/createTimePoints.md),
[`plot.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.TP.md),
[`removeTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/removeTimePoints.md),
[`summary.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.TP.md)

## Examples

``` r
## Create an object of class TP.
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

## Extract the time points from the object.
head(getTimePoints(phenoTP))
#>   timeNumber           timePoint
#> 1          1 2018-05-31 16:37:00
#> 2          2 2018-06-01 09:07:00
#> 3          3 2018-06-01 11:37:00
#> 4          4 2018-06-01 14:37:00
#> 5          5 2018-06-01 16:37:00
#> 6          6 2018-06-02 09:07:00
```
