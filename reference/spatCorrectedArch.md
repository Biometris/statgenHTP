# Maize data corrected for spatial trends.

This dataset contains the corrected data obtained by (1) removing
outliers for single observations and (2) running a spatial model on the
PhenoarchDat1 dataset. See the vignettes for details.

## Usage

``` r
spatCorrectedArch
```

## Format

A data.frame with 37,038 rows and 9 columns:

- timeNumber:

  Time number obtained after formatting the original dataset with the
  function createTP.

- timePoint:

  Original time point.

- LeafArea_corr:

  Leaf area, corrected data

- LeafArea:

  Leaf area from the picture, raw data

- wt:

  Weight factor

- genotype:

  Genotypes

- geno.decomp:

  Combination of treatment levels to decompose the genotypic variance
  (see vignettes)

- colId:

  Column coordinate

- rowId:

  Row coordinate

- plotId:

  Unique pot ID using rowcol coordinates
