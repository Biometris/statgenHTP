#' @export
xyFacetPlot <- function(baseDat,
                        overlayDat = NULL,
                        xVal = "timePoint",
                        yVal,
                        groupVal = "plotId",
                        colVal = groupVal,
                        facetVal = "genotype",
                        title,
                        xLab = "Time",
                        yLab = "Trait") {
  p <- ggplot2::ggplot(plotDat,
                       ggplot2::aes_string(x = xVal, y = yVal,
                                           group = groupVal, color = colVal)) +
    ggplot2::geom_line(show.legend = FALSE) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0, "cm"),
                   panel.border = ggplot2::element_rect(color = "black",
                                                        fill = "transparent"),
                   strip.background = ggplot2::element_rect(color = "black",
                                                            fill = "bisque"),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(title = title, x = xLab, y = yLab)
  if (!is.null(overlayDat)) {
    p + ggplot::geom_line(overlayDat, ggplot2::aes_string(x = xVal, y = yVal),
                          color = "black", size = 2, show.legend = FALSE)
  }
  for (i in 1:ceiling(nlevels(plotDat[[facetVal]]) / 25)) {
    plot(p + ggforce::facet_wrap_paginate(facets = facetVal, nrow = 5, ncol = 5,
                                          page = i))
  }
}
