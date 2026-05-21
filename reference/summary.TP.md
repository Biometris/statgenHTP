# Summary function for TP objects

Function for creating a short summary of the contents of a `TP` object.
The summary consists of the name of the experiment, the number of time
points, the first and last time point and the genotypes defined as
checks.

## Usage

``` r
# S3 method for class 'TP'
summary(object, ...)
```

## Arguments

- object:

  An object of class TP.

- ...:

  Ignored.

## Value

No return value, a summary is printed.

## See also

Other functions for data preparation:
[`as.data.frame.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/as.data.frame.TP.md),
[`createTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/createTimePoints.md),
[`getTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/getTimePoints.md),
[`plot.TP()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.TP.md),
[`removeTimePoints()`](https://biometris.github.io/statgenHTP/index.html/reference/removeTimePoints.md)

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
## Create a summary.
summary(phenoTP)
#> phenoTP contains data for experiment Phenovator.
#> 
#> It contains 73 time points.
#> First time point: 2018-05-31 16:37:00 
#> Last time point: 2018-06-18 16:37:00 
#> 
#> The following genotypes are defined as check genotypes: check1, check2, check3, check4.
```
