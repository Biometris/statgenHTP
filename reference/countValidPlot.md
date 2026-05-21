# Count valid observations per plotId for a given trait

Count valid observations per plotId for a given trait.

## Usage

``` r
countValidPlot(TP, trait, plotIds = NULL)
```

## Arguments

- TP:

  An object of class TP.

- trait:

  A character string indicating the trait for which valid observations
  should be counted.

- plotIds:

  A character vector indicating the plotIds for which valid observations
  should be checked. If `NULL` valid observations are counted for all
  plotIds in TP.

## Value

A named numerical vector with he number of valid observations per
plotId.

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
## Count valid observations for EffpsII for a subset of plots.
countValidPlot(phenoTP,
               trait = "EffpsII",
               plotIds = c("c12r22", "c24r41", "c14r32"))
#> c12r22 c24r41 c14r32 
#>     53     16     73 
```
