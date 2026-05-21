# Maize data, genotypic predictions.

This dataset contains the genotypic predictions obtained by (1) removing
outliers for single observations and (2) running a spatial model on the
PhenoarchDat1 dataset. See the vignettes for details.

## Usage

``` r
spatPredArch
```

## Format

A data.frame with 6,120 rows and 6 columns:

- timeNumber:

  Time number obtained after formatting the original dataset with the
  function createTP.

- timePoint:

  Original time point.

- geno.decomp:

  Combination of treatment levels to decompose the genotypic variance
  (see vignettes)

- genotype:

  Genotypes

- predicted.values:

  Biomass, predicted values

- standard.errors:

  Standard errors associated with the prediction
