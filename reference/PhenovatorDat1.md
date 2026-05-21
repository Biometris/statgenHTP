# Growth chamber data for an Arabidopsis experiment in the Phenovator platform.

A dataset containing data from a growth chamber experiment with
Arabidopsis in the Phenovator platform (WUR, Netherlands, Flood et al.
2016). It consists of one experiment with 1,440 plants grown in a growth
chamber. The number of tested genotypes is 192 with 6 to 7 replicates
per genotype. Four reference genotypes were also tested with 15 or 30
replicates. The studied trait is the photosystem II efficiency (EffpsII)
extracted from the pictures over time (van Rooijen et al. 2017). This
dataset was kindly provided by René Boesten and Mark Aarts.

## Usage

``` r
PhenovatorDat1
```

## Format

A data.frame with 103,839 rows and 10 columns:

- Genotype:

  Genotypes

- Basin:

  Table of experiment

- Replicate:

  Block define after sowing for post-blocking. They are not
  full-resolvable blocks.

- Image_pos:

  Position of the camera

- x:

  Row coordinate

- y:

  Column coordinate

- Sowing_Position:

  Unique pot ID

- timepoints:

  time of picture

- EffpsII:

  Efficiency of the photosystem II

- pos:

  Unique pot ID using rowcol coordinates

## References

Rooijen, Roxanne van, Willem Kruijer, René Boesten, Fred A. van Eeuwijk,
Jeremy Harbinson, and Mark G. M. Aarts. 2017. “Natural Variation of
YELLOW SEEDLING1 Affects Photosynthetic Acclimation of Arabidopsis
Thaliana.” Nature Communications 8 (1).
[doi:10.1038/s41467-017-01576-3](https://doi.org/10.1038/s41467-017-01576-3)

Flood, Pádraic J., Willem Kruijer, Sabine K. Schnabel, Rob van der
Schoor, Henk Jalink, Jan F. H. Snel, Jeremy Harbinson, and Mark G. M.
Aarts. 2016. “Phenomics for Photosynthesis, Growth and Reflectance in
Arabidopsis Thaliana Reveals Circadian and Long-Term Fluctuations in
Heritability.” Plant Methods 12 (1): 14.
[doi:10.1186/s13007-016-0113-y](https://doi.org/10.1186/s13007-016-0113-y)
