# Summary function for fitMod objects

Function for creating a short summary of the contents of a TP object.
The summary consists of the name of the experiment, the number of time
points, the engine used to fit the models and, in case spatial models
where fitted using asreml, the selected spatial model.

## Usage

``` r
# S3 method for class 'fitMod'
summary(object, ...)
```

## Arguments

- object:

  An object of class fitMod.

- ...:

  Ignored.

## Value

No return value, a summary is printed.

## See also

Other functions for spatial modeling:
[`fitModels()`](https://biometris.github.io/statgenHTP/index.html/reference/fitModels.md),
[`getCorrected()`](https://biometris.github.io/statgenHTP/index.html/reference/getCorrected.md),
[`getEffDims()`](https://biometris.github.io/statgenHTP/index.html/reference/getEffDims.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getHerit()`](https://biometris.github.io/statgenHTP/index.html/reference/getHerit.md),
[`getVar()`](https://biometris.github.io/statgenHTP/index.html/reference/getVar.md),
[`plot.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.fitMod.md)

## Examples

``` r
# \donttest{
## Using the first example dataset (PhenovatorDat1):
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

## Fit a SpATS model on few time points:
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = c(1, 6, 36))
#> 2018-05-31 16:37:00
#> 2018-06-02 09:07:00
#> 2018-06-09 14:37:00

## Create a summary.
summary(modPhenoSp)
#> Models in modPhenoSp where fitted for experiment Phenovator.
#> 
#> It contains 3 time points.
#> The models were fitted using SpATS.
#> 
# }
```
