.onLoad <- function(libname = find.package("statgenHTP"),
                    pkgname = "statgenHTP") {
  # CRAN Note avoidance
  if (getRversion() >= "2.15.1")
    utils::globalVariables("i")
  invisible()
}
