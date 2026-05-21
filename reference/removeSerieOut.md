# Replace outliers for series of observations by NA

Function for replacing outliers for series of observations in the data
by NA. The input can either be a data.frame, specified in `dat`, or the
output of the `fitSpline` function, specified in `fitSpline`. Exactly
one of these should be provided as input for the function.

## Usage

``` r
removeSerieOut(
  dat = NULL,
  fitSpline = NULL,
  serieOut,
  reason = c("mean corr", "angle", "slope"),
  traits = attr(x = serieOut, which = "trait")
)
```

## Arguments

- dat:

  A `data.frame`.

- fitSpline:

  An object of class `HTPSpline`, the output of the
  [`fitSpline`](https://biometris.github.io/statgenHTP/index.html/reference/fitSpline.md)
  function.

- serieOut:

  A data.frame with at least the column plotId with values corresponding
  to those in dat/fitSpline.

- reason:

  A character vector indicating which types of outliers should be
  replaced by NA.

- traits:

  The traits that should be replaced by NA. When using the output of
  `detectSerieOut` as input for `serieOut` this defaults to the trait
  used for when detecting the outliers.

## Value

Depending on the input either a `data.frame` or an object of class
`HTPSpline` for which the outliers specified in `serieOut` are replaced
by NA.

## See also

Other functions for detecting outliers for series of observations:
[`detectSerieOut()`](https://biometris.github.io/statgenHTP/index.html/reference/detectSerieOut.md),
[`plot.serieOut()`](https://biometris.github.io/statgenHTP/index.html/reference/plot.serieOut.md)

## Examples

``` r
## Run the function to fit P-splines on a subset of genotypes.
subGenoVator <- c("G160", "G151")
fit.spline <- fitSpline(inDat = spatCorrectedVator,
                        trait = "EffpsII_corr",
                        genotypes = subGenoVator,
                        knots = 50)

## Extract the tables of predicted values and P-spline coefficients.
predDat <- fit.spline$predDat
coefDat <- fit.spline$coefDat

## The coefficients are then used to tag suspect time courses
outVator <- detectSerieOut(corrDat = spatCorrectedVator,
                           predDat = predDat,
                           coefDat = coefDat,
                           trait = "EffpsII_corr",
                           genotypes = subGenoVator,
                           thrCor = 0.9,
                           thrPca = 30,
                           thrSlope = 0.7)

## Replace the outliers by NA in the corrected data.
spatCorrectedVatorOut <- removeSerieOut(dat = spatCorrectedVator,
                                        serieOut = outVator)

## Only replace the slope outliers by NA in the corrected data.
spatCorrectedVatorOut2 <- removeSerieOut(dat = spatCorrectedVator,
                                        serieOut = outVator,
                                        reason = "slope")

## Replace the outliers by NA in the corrected data.
## Replace both the corrected value and the raw trait value by NA.
spatCorrectedVatorOut3 <-
  removeSerieOut(dat = spatCorrectedVator,
                 serieOut = outVator,
                 traits = c("EffpsII", "EffpsII_corr"))
```
