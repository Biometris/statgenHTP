# Extract effective dimensions

Extract effective dimensions from an object of class fitMod. The table
below gives an overview of the effective dimensions and an explanation
of their meaning.

|                     |                                                       |
|---------------------|-------------------------------------------------------|
| Effective Dimension | Explanation                                           |
| colId               | Linear trend along columns                            |
| rowId               | Linear trend along rows                               |
| fCol                | Smooth trend along columns                            |
| fRow                | Smooth trend along rows                               |
| fColRow             | Linear trend in rows changing smoothly along cols     |
| colfRow             | Linear trend in cols changing smoothly along rows     |
| fColfRow            | Smooth-by-smooth interaction trend over rows and cols |
| surface             | Sum of smooth trends                                  |

## Usage

``` r
getEffDims(
  fitMod,
  timePoints = names(fitMod),
  EDType = c("dimension", "ratio"),
  outFile = NULL
)
```

## Arguments

- fitMod:

  An object of class `fitMod`.

- timePoints:

  A character or numeric vector indicating the time point(s) for which
  the effective dimension should be extracted. When using a character
  string to reference a time point, the value has to be an exact match
  to one of the existing time points. When using a number it will be
  matched by its number ("timeNumber") in the timePoints attribute of
  the TP object.

- EDType:

  A character string specifying if the effective dimension ("dimension")
  or the ratio of effective dimensions ("ratio") should be returned.

- outFile:

  A character string indicating the .csv file to which the results
  should be written. If `NULL` no file is written.

## Value

A data.frame with effective dimensions per time point.

## See also

Other functions for spatial modeling:
[`fitModels()`](https://biometris.github.io/statgenHTP/index.html/reference/fitModels.md),
[`getCorrected()`](https://biometris.github.io/statgenHTP/index.html/reference/getCorrected.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getHerit()`](https://biometris.github.io/statgenHTP/index.html/reference/getHerit.md),
[`getVar()`](https://biometris.github.io/statgenHTP/index.html/reference/getVar.md),
[`plot.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.fitMod.md),
[`summary.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.fitMod.md)

## Examples

``` r
# \donttest{
## Using the first example dataset (PhenovatorDat1):
data("PhenovatorDat1")
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

## Fit a SpATS model on few time points:
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = c(1, 6, 20))
#> 2018-05-31 16:37:00
#> 2018-06-02 09:07:00
#> 2018-06-05 14:37:00

## Extract the effective dimensions for all available time points in the
## model object:
effDimSp <- getEffDims(modPhenoSp)
# }
```
