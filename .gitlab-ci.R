pkgsUpdate <- function(quiet = TRUE,
                       instPkgdown = FALSE) {
  install.packages(pkg = "remotes", quiet = quiet)
  if (instPkgdown) {
    install.packages(pkg = "pkgdown", quiet = quiet)
  }
  remotes::install_deps(dependencies = TRUE, quiet = quiet)
  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])
  return(invisible(TRUE))
}


