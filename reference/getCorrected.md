# Extract corrected phenotypic values

Extract corrected phenotypic values from an object of class fitMod.
After fitting a spatial model at each time point, the raw phenotypic
data is corrected by subtracting the (estimated) sources of variation
(environmental, design effect) that are of no interest (nuisances). This
allows keeping the data resolution at the plot/plant level.

## Usage

``` r
getCorrected(fitMod, timePoints = names(fitMod), outFile = NULL)
```

## Arguments

- fitMod:

  An object of class `fitMod`.

- timePoints:

  A character or numeric vector indicating the time point(s) for which
  the corrected values should be extracted. When using a character
  string to reference a time point, the value has to be an exact match
  to one of the existing time points. When using a number it will be
  matched by its number ("timeNumber") in the timePoints attribute of
  the TP object.

- outFile:

  A character string indicating the .csv file to which the results
  should be written. If `NULL` no file is written.

## Value

A data.frame with spatially corrected values per time point.

## See also

Other functions for spatial modeling:
[`fitModels()`](https://biometris.github.io/statgenHTP/index.html/reference/fitModels.md),
[`getEffDims()`](https://biometris.github.io/statgenHTP/index.html/reference/getEffDims.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getHerit()`](https://biometris.github.io/statgenHTP/index.html/reference/getHerit.md),
[`getVar()`](https://biometris.github.io/statgenHTP/index.html/reference/getVar.md),
[`plot.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.fitMod.md),
[`summary.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.fitMod.md)

## Examples

``` r
# \donttest{
## Using the first example dataset (PhenovatorDat1).
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

## Extract the corrected values for one time point:
spatCorrSp <- getCorrected(modPhenoSp,
                           timePoints = 6)
head(spatCorrSp)
#>   timeNumber           timePoint EffpsII_corr EffpsII       wt genotype rowId
#> 1          6 2018-06-02 09:07:00    0.6810552   0.678 1636.287     G001    32
#> 2          6 2018-06-02 09:07:00    0.6625282   0.659 1636.287     G001    58
#> 3          6 2018-06-02 09:07:00    0.6286591   0.630 1636.287     G001    21
#> 4          6 2018-06-02 09:07:00    0.7115547   0.722 1636.287     G001    22
#> 5          6 2018-06-02 09:07:00    0.6764931   0.682 1636.287     G001    33
#> 6          6 2018-06-02 09:07:00    0.6455498   0.641 1636.287     G001     8
#>   colId plotId
#> 1    14 c14r32
#> 2    17 c17r58
#> 3    20 c20r21
#> 4     6  c6r22
#> 5     5  c5r33
#> 6    21  c21r8
# }
```
