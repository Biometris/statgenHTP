# Fit spatial models per time point

Perform REML analysis at each time point using either SpATS or asreml.
The idea is to is to accurately separate the genetic effects from the
spatial effects at each time point. SpATS is used as a default method.
See details for the exact models fitted.

## Usage

``` r
fitModels(
  TP,
  trait,
  timePoints = names(TP),
  extraFixedFactors = NULL,
  geno.decomp = NULL,
  what = c("random", "fixed"),
  useCheck = FALSE,
  useRepId = FALSE,
  engine = c("SpATS", "asreml"),
  spatial = FALSE,
  quiet = FALSE
)
```

## Arguments

- TP:

  An object of class `TP`.

- trait:

  A character string indicating the trait used as response variable in
  the model.

- timePoints:

  A character or numeric vector indicating the time points to be
  modeled. When using a character string to reference a time point, the
  value has to be an exact match to one of the existing time points.
  When using a number it will be matched by its number ("timeNumber") in
  the timePoints attribute of the `TP` object.

- extraFixedFactors:

  A character vector indicating the variables to use as extra fixed
  effects in the model.

- geno.decomp:

  A character vector indicating the variables to use to group the
  genotypic variance in the model.

- what:

  A character vector specifying whether "genotype" should be fitted as
  "random" or "fixed" effect. Note that when using `geno.decomp`,
  fitting a model with genotype as "fixed" effect is not possible.

- useCheck:

  Should check genotypes be used as an extra factor in the model?

- useRepId:

  Should repId be used as a fixed effect in the model? When fitting a
  spatial model rowId and colId are also nested within repId in the
  random part of the model.

- engine:

  A character string indicating the engine used to fit the models.

- spatial:

  Should a spatial model be fitted for asreml?

- quiet:

  Should printed progress messages be suppressed?

## Value

An object of class `fitMod`, a list of fitted models.

## Details

The actual model fitted depends on the function parameters specified.
The basic model is the following:  
trait = **genotype** + e  
In case `useCheck = TRUE`, instead of genotype, genoCheck is used as
genotype and check is used as an extra fixed effect. So then the model
becomes:  
trait = *check* + **genoCheck** + e  
Variables in `extraFixedFactors` are fitted as extra fixed effects.  
  
When `SpATS` is used for modeling, an extra spatial term is always
included in the model. This term is constructed using the function
[`PSANOVA`](https://rdrr.io/pkg/SpATS/man/PSANOVA.html) from the SpATS
package as  
`PSANOVA(colNum, rowNum, nseg = nSeg, nest.div = 2)` where  
`nSeg = c(number of columns, number of rows)`.  
  
When `asreml` is used for modeling and `spatial = TRUE`, four models are
fitted with different covariance structures. The best model is
determined based on a goodness-of-fit criterion, AIC, on 20% of the time
points or at least 10 time points. The best model is then run on all
time points. The following combinations of random and spatial terms are
fitted

- random = repId:rowId + repId:colId, spatial = NULL

- random = repId:rowId + repId:colId, spatial = ar1(rowId):colId

- random = repId:colId + repId:colId, spatial = rowId:ar1(colId)

- random = repId:rowId + repId:colId, spatial = ar1(rowId):ar1(colId)

If there are no replicates in the model, repId is left out from the
random parts above.  
  
When `geno.decomp` is specified, the genotypic variance is decomposed
following the variable(s) chosen. For example, when a treatment is used
in `geno.decomp`, the initial model becomes:  
trait = *treatment* + **treatment:genotype** + e  

## References

Maria Xose Rodriguez-Alvarez, Martin P. Boer, Fred A. van Eeuwijk, Paul
H.C. Eilers (2017). Correcting for spatial heterogeneity in plant
breeding experiments with P-splines. Spatial Statistics
[doi:10.1016/j.spasta.2017.10.003](https://doi.org/10.1016/j.spasta.2017.10.003)
Butler, D. G., et al. (2018). ASReml-R Reference Manual Version 4. VSN
International Ltd, http://asreml.org

## See also

Other functions for spatial modeling:
[`getCorrected()`](https://biometris.github.io/statgenHTP/index.html/reference/getCorrected.md),
[`getEffDims()`](https://biometris.github.io/statgenHTP/index.html/reference/getEffDims.md),
[`getGenoPred()`](https://biometris.github.io/statgenHTP/index.html/reference/getGenoPred.md),
[`getHerit()`](https://biometris.github.io/statgenHTP/index.html/reference/getHerit.md),
[`getVar()`](https://biometris.github.io/statgenHTP/index.html/reference/getVar.md),
[`plot.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.fitMod.md),
[`summary.fitMod()`](https://biometris.github.io/statgenHTP/index.html/reference/summary.fitMod.md)

## Examples

``` r
## Using the first example dataset (PhenovatorDat1):
## Fit a model using SpATS on few time points:
# \donttest{
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
## Fit a model with SpATS for three time points.
modPhenoSp <- fitModels(TP = phenoTP,
                        trait = "EffpsII",
                        timePoints = c(3, 6, 20))
#> 2018-06-01 11:37:00
#> 2018-06-02 09:07:00
#> 2018-06-05 14:37:00
summary(modPhenoSp)
#> Models in modPhenoSp where fitted for experiment Phenovator.
#> 
#> It contains 3 time points.
#> The models were fitted using SpATS.
#> 

## Fit a model with SpATS for a single time point with extra fixed factors
## and check genotypes:
modPhenoSpCheck <- fitModels(TP = phenoTP,
                             trait = "EffpsII",
                             extraFixedFactors = c("repId", "Image_pos"),
                             useCheck = TRUE,
                             timePoints = 3)
#> 2018-06-01 11:37:00


## Fit a model with asreml on few time points with a spatial function:
if (requireNamespace("asreml", quietly = TRUE)) {
  modPhenoSpAs <- fitModels(TP = phenoTP,
                            trait = "EffpsII",
                            timePoints = c(1, 6, 20),
                            engine = "asreml",
                            spatial = TRUE)
}

## Using the second example dataset (PhenoarchDat1):
## Fit a model with SpATS on one time points with two variables for
## geno.decomp:
data("PhenoarchDat1")
phenoTParch <- createTimePoints(dat = PhenoarchDat1,
                                experimentName = "Phenoarch",
                                genotype = "Genotype",
                                timePoint = "Date",
                                plotId = "pos",
                                rowNum = "Row",
                                colNum = "Col")

modPhenoSpGD <- fitModels(TP = phenoTParch,
                          trait = "LeafArea",
                          geno.decomp = c("Scenario", "population"),
                          timePoints = 16)
#> 2017-04-28
# }

```
