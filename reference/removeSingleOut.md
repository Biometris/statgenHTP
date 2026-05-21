# Replace outliers for single observations by NA

Function for replacing outliers for single observations by NA.

## Usage

``` r
removeSingleOut(TP, singleOut, trait = attr(x = singleOut, which = "trait"))
```

## Arguments

- TP:

  An object of class TP.

- singleOut:

  A data.frame with at least the columns plotId and timePoint with
  values corresponding to those in TP. If a column outlier is present,
  as in the output of `detectSingleOut`, only plotId x timePoint
  combinations for which outlier = 1 will be set to NA. If no column
  outlier is present, all observations in singleOut will be set to NA.

- trait:

  The trait that should be set to NA. Can be ignored when using the
  output of `detectSingleOut` as input.

## Value

An object of class TP, the input with the outlier replaced by NA.

## See also

Other functions for detecting outliers for single observations:
[`detectSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOut.md),
[`detectSingleOutMaize()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOutMaize.md),
[`plot.singleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.singleOut.md)

## Examples

``` r
## Create a TP object containing the data from the Phenovator.
PhenovatorDat1 <- PhenovatorDat1[!PhenovatorDat1$pos %in%
                                 c("c24r41", "c7r18", "c7r49"), ]
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

## First select a subset of plants, for example here 9 plants.
plantSel <- phenoTP[[1]]$plotId[1:9]
# Then run on the subset
resuVatorHTP <- detectSingleOut(TP = phenoTP,
                                trait = "EffpsII",
                                plotIds = plantSel,
                                confIntSize = 3,
                                nnLocfit = 0.1)

## Replace the studied trait by NA for the plants marked as outliers.
phenoTPOut <- removeSingleOut(phenoTP, resuVatorHTP)
```
