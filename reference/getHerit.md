# Extract heritabilities

Extract heritabilities from an object of class fitMod. When
`geno.decomp` is used, the heritabilities of each level of geno.decomp
are stored in separate columns.

## Usage

``` r
getHerit(fitMod, timePoints = names(fitMod), outFile = NULL)
```

## Arguments

- fitMod:

  An object of class `fitMod`.

- timePoints:

  A character or numeric vector indicating the time point(s) for which
  the heritabilities should be extracted. When using a character string
  to reference a time point, the value has to be an exact match to one
  of the existing time points. When using a number it will be matched by
  its number ("timeNumber") in the timePoints attribute of the TP
  object.

- outFile:

  A character string indicating the .csv file to which the results
  should be written. If `NULL` no file is written.

## Value

A data.frame with heritabilities per time point.

## See also

Other functions for spatial modeling:
[`fitModels()`](https://biometris.github.io/statgenHTP/index.html/reference/fitModels.md),
[`getCorrected()`](https://biometris.github.io/statgenHTP/index.html/reference/getCorrected.md),
[`getEffDims()`](https://biometris.github.io/statgenHTP/index.html/reference/getEffDims.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getVar()`](https://biometris.github.io/statgenHTP/index.html/reference/getVar.md),
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

## Extract the heritabilities for all available time points.#'
getHerit(modPhenoSp)
#>   timeNumber           timePoint   h2
#> 1          1 2018-05-31 16:37:00 0.70
#> 2          6 2018-06-02 09:07:00 0.71
#> 3         20 2018-06-05 14:37:00 0.88
# }
```
