# Arabidopsis data corrected for spatial trends.

This dataset contains the corrected data obtained by (1) removing
outliers for single observations and (2) running a spatial model on the
PhenovatorDat1 dataset. See the vignettes for details.

## Usage

``` r
spatCorrectedVator
```

## Format

A data.frame with 103,801 rows and 11 columns:

- timeNumber:

  Time number obtained after formatting the original dataset with the
  function createTP.

- timePoint:

  Original time point.

- EffpsII_corr:

  Efficiency of the photosystem II, corrected data

- EffpsII:

  Efficiency of the photosystem II, raw data

- genotype:

  Genotypes

- repId:

  Block define after sowing for post-blocking.

- Image_pos:

  Position of the camera

- check:

  Status of the genotypes: check for the reference genotypes, noCheck
  for the others.

- colId:

  Column coordinate

- rowId:

  Row coordinate

- plotId:

  Unique pot ID using rowcol coordinates
