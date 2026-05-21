# Extract variances

Extract variances from an object of class fitMod.

## Usage

``` r
getVar(fitMod, timePoints = names(fitMod), outFile = NULL)
```

## Arguments

- fitMod:

  An object of class `fitMod`.

- timePoints:

  A character or numeric vector indicating the time point(s) for which
  the variances should be extracted. When using a character string to
  reference a time point, the value has to be an exact match to one of
  the existing time points. When using a number it will be matched by
  its number ("timeNumber") in the timePoints attribute of the TP
  object.

- outFile:

  A character string indicating the .csv file to which the results
  should be written. If `NULL` no file is written.

## Value

A data.frame with variances per time point.

## See also

Other functions for spatial modeling:
[`fitModels()`](https://biometris.github.io/statgenHTP/index.html/reference/fitModels.md),
[`getCorrected()`](https://biometris.github.io/statgenHTP/index.html/reference/getCorrected.md),
[`getEffDims()`](https://biometris.github.io/statgenHTP/index.html/reference/getEffDims.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getHerit()`](https://biometris.github.io/statgenHTP/index.html/reference/getHerit.md),
[`plot.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.fitMod.md),
[`summary.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.fitMod.md)

## Examples

``` r
# \donttest{
## Using the first example dataset (PhenovatorDat1):
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

## Fit a SpATS model on few time points.
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = c(1, 6, 20))
#> 2018-05-31 16:37:00
#> 2018-06-02 09:07:00
#> 2018-06-05 14:37:00

## Extract the variances for all available time points.
getVar(modPhenoSp)
#>   timeNumber           timePoint       varGen       varRes       varCol
#> 1          1 2018-05-31 16:37:00 0.0002565400 0.0007483918 1.199643e-05
#> 2          6 2018-06-02 09:07:00 0.0002019810 0.0005521569 6.760718e-07
#> 3         20 2018-06-05 14:37:00 0.0002663746 0.0002319816 5.434449e-06
#>         varRow
#> 1 8.493845e-05
#> 2 5.966676e-05
#> 3 6.425597e-05
# }
```
