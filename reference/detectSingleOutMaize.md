# detectSingleOutMaize

Function to detect plant outliers in a temporal lattice experiment on
Maize which can be extended to other experiment types. The criteria
needs three phenotypes (ex for maize: the estimated biomass, plant
height and phyllocron)

- plants are identified as "small outlier plant":

  if for biomass AND phyllocron \\res_i \< \mu\_{res} - qnorm(threshold)
  \* sd\_{res}\\

- plants are identified as "big outlier plant":

  if for biomass AND plant height \\res_i \> \mu\_{res} +
  qnorm(threshold) \* sd\_{res}\\

## Usage

``` r
detectSingleOutMaize(
  TP,
  timeBeforeTrt,
  trait1 = "Biomass",
  trait2 = "PlantHeight",
  trait3 = "phyllocron",
  thr = 0.95
)
```

## Arguments

- TP:

  An object of class TP.

- timeBeforeTrt:

  A character or numeric value indicating the date just before treatment
  in the experiment. When using a character string to reference a time
  point, the value has to be an exact match to one of the existing
  timePoints. When using a number it will be matched by its number
  ("timeNumber") in the timePoints attribute of the TP object.

- trait1:

  A character vector indicating the first trait to model in TP.

- trait2:

  A character vector indicating the second trait to model in TP.

- trait3:

  A character vector indicating the third trait to model in TP.

- thr:

  A numeric value indicating the threshold.

## Value

A list with three data.frames, `modDat` containing the data used for
fitting the models, `smallPlants` containing the plants identified as
small plants and `bigPlants` containing the plants identified as big
plants.

## See also

Other functions for detecting outliers for single observations:
[`detectSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSingleOut.md),
[`plot.singleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.singleOut.md),
[`removeSingleOut()`](https://biometris.github.io/statgenHTP/index.html/reference/removeSingleOut.md)

## Examples

``` r
# \donttest{
## Create a TP object containing the data from the PhenoArch.
phenoTParch <- createTimePoints(dat = PhenoarchDat1,
                                experimentName = "Phenoarch",
                                genotype = "Genotype",
                                timePoint = "Date",
                                plotId = "pos",
                                rowNum = "Row",
                                colNum = "Col")
singleOutMaize <- detectSingleOutMaize(phenoTParch,
                                       timeBeforeTrt = "2017-04-27",
                                       trait1 = "Biomass",
                                       trait2 = "PlantHeight",
                                       trait3 = "phyllocron",
                                       thr = 0.95)
# }
```
