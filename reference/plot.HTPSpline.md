# Plot the results of a fitted spline.

Plot the results of a fitted spline.

## Usage

``` r
# S3 method for class 'HTPSpline'
plot(
  x,
  ...,
  plotType = c("predictions", "derivatives", "derivatives2"),
  genotypes = NULL,
  plotIds = NULL,
  title = NULL,
  output = TRUE,
  outFile = NULL,
  outFileOpts = NULL
)
```

## Arguments

- x:

  An object of class `HTPSpline`.

- ...:

  Ignored.

- plotType:

  A character string indicating which spline component should be
  plotted, either predictions, derivatives or second derivatives
  ("derivatives2").

- genotypes:

  A character vector indicating the genotypes for which spline
  components should be plotted.

- plotIds:

  A character vector indicating the plotIds for which spline components
  should be plotted.

- title:

  A character string used as title for the plot. If `NULL` a default
  title is added to the plot depending on `plotType`.

- output:

  Should the plot be output to the current device? If `FALSE` only a
  (list of) ggplot object(s) is invisibly returned. Ignored if `outFile`
  is specified.

- outFile:

  A character string indicating the .pdf file to which the plots should
  be written. If `NULL`, no file is written.

- outFileOpts:

  A named list of extra options for the pdf outfile, e.g. width and
  height. See [`pdf`](https://rdrr.io/r/grDevices/pdf.html) for all
  possible options.

## Value

A list of object of class ggplot is invisibly returned.

## See also

Other functions for fitting splines:
[`fitSpline()`](https://biometris.github.io/statgenHTP/index.html/reference/fitSpline.md)

## Examples

``` r
## The data from the Phenovator platform have been corrected for spatial
## trends and outliers for single observations have been removed.

## Fit P-Splines on a subset of genotypes
subGeno <- c("G070", "G160")
fit.spline <- fitSpline(inDat = spatCorrectedVator,
                        trait = "EffpsII_corr",
                        genotypes = subGeno,
                        knots = 50)

## Visualize the P-Spline predictions for one genotype.
plot(fit.spline, genotypes = "G160")


## Visualize the first and second derivatives of the predictions for one plant.
plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives")

plot(fit.spline, plotIds = "c10r29", plotType =  "derivatives2")

```
