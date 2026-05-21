# Count valid observations per time point for a given trait

Count valid observations per time point for a given trait.

## Usage

``` r
countValid(TP, trait)
```

## Arguments

- TP:

  An object of class TP.

- trait:

  A character string indicating the trait for which valid observations
  should be counted.

## Value

A named numerical vector with he number of valid observations per time
point .

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
## Count valid observations for EffpsII per time point.
validPheno <- countValid(phenoTP, trait = "EffpsII")
head(validPheno)
#> 2018-05-31 16:37:00 2018-06-01 09:07:00 2018-06-01 11:37:00 2018-06-01 14:37:00 
#>                1410                1391                1407                1413 
#> 2018-06-01 16:37:00 2018-06-02 09:07:00 
#>                1413                1411 
```
