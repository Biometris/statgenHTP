#' Helper function for prediction. It returns the design matrix for the B-splines
#' defined by the same knots that are used in the fitting stage.
#'
#' @noRd
#' @keywords internal
spline.bbase <- function(knots,
                         X.,
                         BDEG.,
                         deriv = 0) {
  B <- splines::spline.des(knots = knots, x = X., derivs = deriv,
                           ord = BDEG. + 1, outer.ok = TRUE)$design
  return(B)
}

#' Helper function to compute a B-spline basis matrix using evenly spaced knots.
#'
#' @noRd
#' @keywords internal
bbase <- function(x,
                  xl = min(x),
                  xr = max(x),
                  ndx = 10,
                  bdeg = 3) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::spline.des(knots = knots, x = x, ord = bdeg + 1,
                           outer.ok = TRUE)$design
  res <- list(B = B, knots = knots)
  return(res)
}

#' Helper function to calculate the singular value decomposition of the B-spline
#' basis matrix to obtain the mixed model design matrices.
#'
#' @noRd
#' @keywords internal
MM.basis <- function (x,
                      xl = min(x),
                      xr = max(x),
                      ndx,
                      bdeg,
                      pord) {
  Bb <- bbase(x = x, xl = xl, xr = xr, ndx = ndx, bdeg = bdeg)
  knots <- Bb$knots
  B <- Bb$B
  m <- ncol(B)
  n <- nrow(B)
  D <- diff(diag(m), differences = pord)
  P.svd <- svd(crossprod(D))
  d <- (P.svd$d)[1:(m - pord)]
  U.Z <- (P.svd$u)[, 1:(m - pord)]
  Z <- Matrix::Matrix(B %*% U.Z)
  X <- Matrix::Matrix(outer(X = x, Y = 0:(pord-1), FUN = "^"))
  U.X <- outer(X = knots[-c((1:(bdeg - 1)),
                            (length(knots) - bdeg + 2):length(knots))],
               Y = 0:(pord-1), FUN = "^")
  res <- list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots,
              U.X = U.X, U.Z = U.Z)
  return(res)
}

#' Helper function to construct the components of the precision matrix
#' (as needed by the algorithm).
#'
#' @noRd
#' @keywords internal
constructCapitalLambda <- function(g) {
  length.eq <- all(sapply(X = g, FUN = function(x) {
    diff(range(unlist(lapply(X = x, FUN = length)))) < .Machine$double.eps ^ 0.5
  }))
  if (length.eq) {
    l <- length(g)
    if (l == 1) {
      if (length(g[[1]]) == 1) {
        res <- g
      } else {
        res <- do.call("c", lapply(X = g, FUN = function(x) x))
      }
    } else {
      dim <- sapply(X = g, FUN = function(x) {
        if (is.list(x)) {
          unlist(lapply(X = x, FUN = length))[1]
        } else {
          length(x)
        }
      })
      end <- cumsum(dim)
      init <- end - dim + 1
      res <- do.call("c", lapply(X = 1:length(g),
                                 FUN = function(x, g, init, end, dim) {
                                   temp <- g[[x]]
                                   if (is.list(temp)) {
                                     lapply(X = temp, FUN = function(y, x, dim) {
                                       aux <- rep(0L, length.out = sum(dim))
                                       aux[init[x]:end[x]] <- y
                                       return(aux)
                                     }, x = x, dim = dim)
                                   } else {
                                     aux <- rep(0L, length.out = sum(dim))
                                     aux[init[x]:end[x]] <- temp
                                     list(aux)
                                   }
                                 }, g = g, init = init, end = end, dim = dim))
    }
  } else {
    stop("Error in constructCapitalLambda")
  }
  return(res)
}

#' Helper function to construct the Henderson's equations.
#'
#' @noRd
#' @keywords internal
construct.block <- function(A1,
                            A2,
                            A3,
                            A4) {
  block <- rbind(cbind(A1, A2), cbind(A3, A4))
  return(block)
}

#' Helper function to compute the row-wise Kronecker product of one matrix as
#' indicated in GLAM (Eilers et al., 2006).
#'
#' @noRd
#' @keywords internal
Rten <- function(X) {
  one <- Matrix::Matrix(1, 1, ncol(X))
  kronecker(X, one) * kronecker(one, X)
}

#' Helper function to compute the row-wise Kronecker product of two matrices as
#' indicated in GLAM (Eilers et al., 2006).
#'
#' @noRd
#' @keywords internal
Rten2 <- function(X1,
                  X2) {
  one.1 <- Matrix::Matrix(1, 1, ncol(X1))
  one.2 <- Matrix::Matrix(1, 1, ncol(X2))
  kronecker(X1, one.2) * kronecker(one.1, X2)
}

#' Helper function for array dimension rotation as indicated in GLAM
#' (Eilers et al., 2006).
#'
#' @noRd
#' @keywords internal
RH <- function(X,
               A) {
  Matrix::t(X %*% A)
}

#' Helper function for fast computation of the elements of the
#' principal diagonal of the Henderson system of equations using GLAM
#' (Eilers et al., 2006)
#'
#' @importFrom methods as new
#'
#' @noRd
#' @keywords internal
A1.form <- function(l,
                    w = NULL) {
  d <- length(l)
  n <- rev(sapply(X = l, FUN = nrow))
  c1 <- rev(sapply(X = l, FUN = ncol))
  if (is.null(w)) {
    W <- Matrix::Matrix(1, nrow = n[1], ncol = n[2])
  } else {
    dim(w) <- n
    W <- as(w, "dgeMatrix")
  }
  lRTen <- lapply(X = rev(l), FUN = Rten)
  tmp <- Reduce(Matrix::crossprod, x = lRTen, init = W)
  tmp <- array(tmp, dim = rep(c1, rep(2, d)))
  Fast <- if (prod(c1)) {
    bdiag_m(lapply(X = seq_len(dim(tmp)[3]), FUN = function(i) {
      tmp[, , i, i]
    }))
  } else {
    aperm(tmp, c(2 * (1:d) - 1, 2 * (1:d)))
  }
  return(Fast)
}

#' Helper function for for fast computation of the elements of the
#' secondary diagonal of the Henderson system of equations using GLAM
#' (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
A2.form <- function(l1,
                    l2,
                    w = NULL) {
  d1 <- length(l1)
  d2 <- length(l2)
  if (d1 != d2) {
    stop("l1 and l2 should have the same dimension")
  }
  n <- rev(sapply(X = l1, FUN = nrow))
  d <- rev(sapply(X = l1, FUN = ncol))
  c1 <- rev(sapply(X = l2, FUN = ncol))
  if (is.null(w)) {
    W <- Matrix::Matrix(1, nrow = n[1], ncol = n[2])
  } else {
    dim(w) <- n
    W <- as(w, "dgeMatrix")
  }
  lRTen2 <- mapply(FUN = Rten2, rev(l2), rev(l1))
  tmp <- Reduce(Matrix::crossprod, x = lRTen2, init = W)
  tmp <- array(tmp, dim = as.vector(rbind(d, c1)))
  Fast1 <- aperm(tmp, c(2 * (1:d1) - 1, 2 * (1:d1)))
  Fast <- if (prod(d)) {
    Matrix::Matrix(Fast1, nrow = prod(d), sparse = TRUE)
  } else {
    aperm(tmp, c(2 * (1:d1) - 1, 2 * (1:d1)))
  }
  return(Fast)
}

#' Helper function for fast computation of the XtX element
#' of the Henderson system of equations using GLAM (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
XtX <- function(X,
                w = NULL) {
  A1.form(X, w)
}

#' Helper function for fast computation of the XtZ element
#' of the Henderson system of equations using GLAM (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
XtZ <- function(X,
                Z,
                w = NULL) {
  d <- length(Z)
  res <- NULL
  for (i in 1:d) {
    res <- cbind(res, A2.form(X, Z[[i]], w))
  }
  res
}

#' Helper function for fast computation of the ZtZ element
#' of the Henderson system of equations using GLAM (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
ZtZ <- function(Z,
                w = NULL) {
  d <- length(Z)
  upper <- vector(mode = "list", length = d)
  for (i in 1:d) {
    upper[[i]] <- vector(mode = "list", length = d)
    upper[[i]][[i]] <- A1.form(Z[[i]], w)
  }
  ## Obtain the elements of the matrix.
  if (d > 1) {
    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        upper[[i]][[j]] <- A2.form(Z[[i]], Z[[j]], w)
      }
    }
  }
  ## Create the matrix.
  res <- Reduce(cbind, x = upper[[1]])
  for (i in 2:d) {
    tmp <- Reduce(cbind, x = upper[[i]])
    for (j in (i - 1):1) {
      if (length(upper[[j]][[i]])) {
        tmp <- cbind(Matrix::t(upper[[j]][[i]]), tmp)
      }
    }
    if (nrow(tmp)) {
      res <- rbind(res, tmp)
    }
  }
  return(res)
}

#' Helper function for fast computation of the Xty element
#' of the Henderson system of equations using GLAM (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
Xty <- function(X,
                y,
                w = NULL) {
  nRow <- nrow(utils::tail(X, n = 1)[[1]])
  if (is.null(w)) {
    Y <- Matrix::Matrix(y, nrow = nRow)
  } else {
    Y <- Matrix::Matrix(w * y, nrow = nRow)
  }
  tmp <- Reduce(Matrix::crossprod, x = rev(X), init = Y)
  as.vector(tmp)
}

#' Helper function for fast computation of the Zty element
#' of the Henderson system of equations using GLAM (Eilers et al., 2006)
#'
#' @noRd
#' @keywords internal
Zty <- function(Z,
                y,
                w = NULL) {
  nRow <- nrow(utils::tail(Z[[1]], n = 1)[[1]])
  if (is.null(w)) {
    Y <- Matrix::Matrix(y, nrow = nRow)
  } else {
    Y <- Matrix::Matrix(w * y, nrow = nRow)
  }
  res <- unlist(lapply(X = Z, FUN = function(z) {
    as.vector(Reduce(Matrix::crossprod, x = rev(z), init = Y))
  }), use.names = FALSE)
  return(res)
}

#' Helper function to compute the estimated values for the fixed effects using
#' GLAM (Eilers et al., 2006).
#'
#' @noRd
#' @keywords internal
Xtheta <- function(X,
                   theta) {
  d <- length(X)
  n <- rev(sapply(X = X, FUN = ncol))
  Theta <- Matrix::Matrix(theta, nrow = n[1])
  tmp <- RH(X[[d]], Theta)
  for (i in (d - 1):1) {
    tmp <- RH(X[[i]], tmp)
  }
  as.vector(tmp)
}

#' Helper function to compute the estimated values for the random effects
#' using GLAM (Eilers et al., 2006).
#'
#' @noRd
#' @keywords internal
Ztheta <- function(Z,
                   theta,
                   np) {
  d <- length(Z)
  for (i in 1:d) {
    if (i == 1) {
      res <- Xtheta(Z[[i]], theta[1:(np[1])])
    } else {
      init <- sum(np[1:(i - 1)])
      fin <- np[i]
      if (fin) {
        res <- res + Xtheta(Z[[i]], theta[(init + 1):(init + fin)])
      }
    }
  }
  res
}

#' Helper function to construct the Henderson system of equations using GLAM
#' (Eilers et al., 2006), or traditional matrix arithmetic.
#'
#' @noRd
#' @keywords internal
construct.matrices <- function(X,
                               Z,
                               z,
                               w,
                               GLAM) {
  if (GLAM) {
    XtX. <- XtX(X, w)
    XtZ. <- XtZ(X, Z, w)
    ZtX. <- Matrix::t(XtZ.)
    ZtZ. <- ZtZ(Z, w)
    Zty. <- Zty(Z, z, w)
    Xty. <- Xty(X, z, w)
    yty. <- sum((z ^ 2) * w)
    ZtXtZ <- rbind(XtZ., ZtZ.)
    u <- c(Xty., Zty.)
  } else {
    XtW. <- t(X * w)
    XtX. <- XtW. %*% X
    XtZ. <- XtW. %*% Z
    ZtX. <- Matrix::t(XtZ.)
    ZtW. <- Matrix::t(Z * w)
    ZtZ. <- ZtW. %*% Z
    Xty. <- XtW. %*% z
    Zty. <- ZtW. %*% z
    yty. <- sum((z ^ 2) * w)
    ZtXtZ = rbind(XtZ., ZtZ.)
    u <- c(Xty., Zty.)
  }
  res <- list(XtX. = XtX., XtZ. = XtZ., ZtX. = ZtX., ZtZ. = ZtZ., Xty. = Xty.,
              Zty. = Zty., yty. = yty., ZtXtZ = ZtXtZ, u = u)
  return(res)
}

#' Helper function to compute the REML deviance.
#'
#' @noRd
#' @keywords internal
deviance_spam <- function(C,
                          G,
                          w,
                          sigma2,
                          ssr,
                          edf) {
  log_det_C <- spam::determinant(C)$modulus * 2
  log_det_G <- spam::determinant(G)$modulus
  deviance <- log_det_C + log_det_G + sum(log(sigma2 / w)) + ssr / sigma2 + edf
  return(deviance)
}
