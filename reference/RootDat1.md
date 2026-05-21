# Greenhouse data for an experiment in the RootPhAir platform.

A dataset containing greenhouse data from the RootPhAir platform
(UCLouvain, Belgium). It consists of one experiment with one aeroponic
tanks with 340 maize plants. The studied traits are the root tip
coordinates in y and x axis, extracted from the pictures over time.
Plants were pictured every 2 hours for 10 days. This dataset was kindly
provided by Xavier Draye.

## Usage

``` r
RootDat1
```

## Format

A data.frame with 16,275 rows and 10 columns:

- Exp:

  Experiment number

- thermalTime:

  Thermal time cumulated

- Genotype:

  Genotype

- plantId:

  Unique pot using tank and rowcol coordinate

- Tank:

  Tank A or B

- Strip:

  Number of strip of five plants (i.e. row coordinate)

- Pos:

  Position within th strip (i.e. column coordinate)

- tipPos_x:

  Position of the root tip in x axis

- tipPos_y:

  Position of the root tip in y axis

- Time:

  Time of measurement
