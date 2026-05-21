# Detect outliers for single observations

Detect outlying observations in a time series by modeling each plotId
using a local regression.

## Usage

``` r
detectSingleOut(
  TP,
  trait,
  plotIds = NULL,
  checkEdges = TRUE,
  confIntSize = 5,
  nnLocfit = 0.5
)
```

## Arguments

- TP:

  An object of class `TP`.

- trait:

  A character vector indicating the trait to model in `TP`.

- plotIds:

  A character vector of plotIds for which the outliers should be
  detected. If `NULL`, all plotIds in `TP` are used.

- checkEdges:

  Before fitting the local regression should a check be done if the
  first and last time point for a plot are outlying observations?

- confIntSize:

  A numeric value defining the confidence interval (see Details).

- nnLocfit:

  A numeric value defining the constant component of the smoothing
  parameter nn (see Details).

## Value

An object of class singleOut, a `data.frame` with the following columns.

- plotId:

  plotId

- timePoint:

  time point

- trait:

  modeled trait

- yPred:

  prediction from the local regression

- sd_yPred:

  standard deviation of the prediction

- lwr:

  lower bound of the confidence interval

- upr:

  upper bound of the confidence interval

- outlier:

  flag for detected outlier (a value of 1 indicates the observation is
  an outlier)

## Details

See locfit() help function from the locfit R library. The user can act
on:

- nnLocfit:

  the constant of the smoothing parameter. Increase nnLocfit to have a
  very smooth curve

- confIntSize:

  the level to calculate the confidence interval. Increase confIntSize
  to exclude less outliers

## See also

Other functions for detecting outliers for single observations:
[`detectSingleOutMaize()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOutMaize.md),
[`plot.singleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.singleOut.md),
[`removeSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/removeSingleOut.md)

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

## First select a subset of plants, for example here 9 plants
plantSel <- phenoTP[[1]]$plotId[1:9]
# Then run on the subset
resuVatorHTP <- detectSingleOut(TP = phenoTP,
                                trait = "EffpsII",
                                plotIds = plantSel,
                                confIntSize = 3,
                                nnLocfit = 0.1)
```
