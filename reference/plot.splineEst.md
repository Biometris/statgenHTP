# Plot the results of estimated spline parameters.

Plot the results of estimated spline parameters.

## Usage

``` r
# S3 method for class 'splineEst'
plot(
  x,
  ...,
  plotType = c("box", "hist"),
  what = attr(x, "what"),
  title = NULL,
  output = TRUE,
  outFile = NULL,
  outFileOpts = NULL
)
```

## Arguments

- x:

  An object of class `splineEst`

- ...:

  Ignored.

- plotType:

  A character string indicating the type of plot to be made.

- what:

  The types of estimate that should be plotted.

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

A list of objects of class ggplot is invisibly returned.

## See also

Other functions for spline parameter estimation:
[`estimateSplineParameters()`](https://biometris.github.io/statgenHTP/index.html/reference/estimateSplineParameters.md)
