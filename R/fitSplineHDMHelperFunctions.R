#' Helper function for ....
#'
#' @noRd
#' @keywords internal
spline.bbase <- function(knots,
                         X.,
                         BDEG.,
                         deriv = 0,
                         eps = 1e-15) {
  B <- splines::spline.des(knots = knots, x = X., derivs = deriv,
                           ord = BDEG. + 1, outer.ok = TRUE)$design
  return(B)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
bbase <- function(x,
                  xl = min(x),
                  xr = max(x),
                  ndx = 10,
                  bdeg = 3,
                  eps = 1e-15) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  B <- splines::spline.des(knots = knots, x = x, ord = bdeg + 1,
                           outer.ok = TRUE)$design
  res <- list(B = B, knots = knots)
  return(res)
}

#' Helper function for ....
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
  Z <- B %*% U.Z
  X <- NULL
  for (i in 0:(pord - 1)){
    X <- cbind(X, x ^ i)
  }
  U.X <- NULL
  for (i in 0:(pord - 1)) {
    U.X <- cbind(U.X,
                 knots[-c((1:(bdeg - 1)),
                          (length(knots) - (bdeg - 1) + 1):length(knots))] ^ i)
  }
  res <- list(X = X, Z = Z, d = d, B = B, m = m, D = D, knots = knots,
              U.X = U.X, U.Z = U.Z)
  return(res)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
construct.capital.lambda <- function(g) {
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
    stop("Error in construct.capital.lambda")
  }
  return(res)
}

#' Helper function for ....
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

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Rten <- function(X) {
  one <- Matrix::Matrix(1, 1, ncol(X))
  kronecker(X, one) * kronecker(one, X)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Rten2 <- function(X1,
                  X2) {
  one.1 <- Matrix::Matrix(1, 1, ncol(X1))
  one.2 <- Matrix::Matrix(1, 1, ncol(X2))
  kronecker(X1, one.2) * kronecker(one.1, X2)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
H <- function(X,
              A) {
  d <- dim(A)
  M <- Matrix::Matrix(A, nrow = d[1])
  XM <- X %*% M
  slam::as.simple_sparse_array(array(XM, dim = c(nrow(XM), d[-1])))
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Rotate <- function(A) {
  d <- 1:length(dim(A))
  d1 <- c(d[-1], d[1])
  slam::as.simple_sparse_array(aperm(A, d1))
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
RH <- function(X,
               A) {
  Rotate(H(X, A))
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
A1.form <- function(l,
                    w = NULL) {
  d <- length(l)
  n <- rev(sapply(X = l, FUN = nrow))
  c <- rev(sapply(X = l, FUN = ncol))
  if (is.null(w)) {
    W <- slam::as.simple_sparse_array(array(1, dim = n))
  } else {
    W <- slam::as.simple_sparse_array(array(w, dim = n))
  }
  tmp <- RH(Matrix::t(Rten(l[[d]])), W)
  for (i in (d-1):1) {
    tmp <- RH(Matrix::t(Rten(l[[i]])), tmp)
  }
  dim(tmp)<- rep(c, rep(2,d))
  Fast1 <- slam::as.simple_sparse_array(aperm(tmp, as.vector(
    matrix(1:(d * 2), byrow = TRUE, ncol = 2))))
  Fast <- if (prod(c)) {
    Matrix::Matrix(Fast1, nrow = prod(c))
  } else {
    Fast1
  }
  return(Fast)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
A2.form <- function(l1,
                    l2,
                    w = NULL) {
  d1 <- length(l1)
  d2 <- length(l2)
  if (!(d1 == d2)) {
    stop("l1 and l2 should have the same dimension")
  }
  n <- rev(sapply(X = l1, FUN = nrow))
  d <- rev(sapply(X = l1, FUN = ncol))
  c <- rev(sapply(X = l2, FUN = ncol))
  if (is.null(w)) {
    W <- slam::as.simple_sparse_array(array(1, dim = n))
  } else {
    W <- slam::as.simple_sparse_array(array(w, dim = n))
  }
  tmp <- RH(Matrix::t(Rten2(l2[[d1]], l1[[d1]])), W)
  for (i in (d1-1):1) {
    tmp <- RH(Matrix::t(Rten2(l2[[i]], l1[[i]])), tmp)
  }
  dim(tmp)<- as.vector(rbind(d, c))
  Fast1 <- slam::as.simple_sparse_array(aperm(tmp, as.vector(
    matrix(1:(d1 * 2), byrow = TRUE, ncol = 2))))
  Fast <- if (prod(d)) {
    Matrix::Matrix(Fast1, nrow = prod(d))
  } else {
    Fast1
  }
  return(Fast)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
XtX <- function(X,
                w = NULL) {
  A1.form(X, w)
}

#' Helper function for ....
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

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
ZtZ <- function(Z,
                w = NULL) {
  d <- length(Z)
  upper <- list()
  for (i in 1:d) {
    upper[[i]] <- list()
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
  res <- NULL
  for (i in 1:d) {
    if (i == 1) {
      res <- do.call("cbind", upper[[1]])
    } else {
      tmp <- do.call("cbind", upper[[i]])
      for (j in (i - 1):1) {
        if (length(upper[[j]][[i]])) {
          tmp <- cbind(Matrix::t(upper[[j]][[i]]), tmp)
        }
      }
      if (nrow(tmp)) {
        res <- rbind(res, tmp)
      }
    }
  }
  return(res)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Xty <- function(X,
                y,
                w = NULL) {
  d <- length(X)
  n <- rev(sapply(X = X, FUN = nrow))
  if (is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w * y, n)
  }
  tmp <- RH(Matrix::t(X[[d]]), Y)
  for (i in (d - 1):1) {
    tmp <- RH(Matrix::t(X[[i]]), tmp)
  }
  as.vector(tmp)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Zty <- function(Z,
                y,
                w = NULL) {
  d <- length(Z)
  n <- rev(sapply(X = Z[[1]], FUN = nrow))
  if (is.null(w)) {
    Y <- array(y, n)
  } else {
    Y <- array(w * y, n)
  }
  res <- NULL
  for (i in 1:d) {
    k <- length(Z[[i]])
    tmp <- RH(Matrix::t(Z[[i]][[k]]), Y)
    for (j in (k - 1):1) {
      tmp <- RH(Matrix::t(Z[[i]][[j]]), tmp)
    }
    res <- c(res, as.vector(tmp))
  }
  return(res)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
Xtheta <- function(X,
                   theta) {
  d <- length(X)
  n <- rev(sapply(X = X, FUN = ncol))
  Theta <- array(theta, dim = n)
  tmp <- RH(X[[d]], Theta)
  for (i in (d - 1):1) {
    tmp <- RH(X[[i]], tmp)
  }
  as.vector(tmp)
}

#' Helper function for ....
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

#' Helper function for ....
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

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
construct.G.matrix <- function(g,
                               la) {
  la.g <- la[-1]
  g.comp <- length(g)
  g.tot.comp <- unlist(lapply(X = g, FUN = length))
  np.e.g <- cumsum(g.tot.comp)
  np.s.g <- np.e.g - g.tot.comp + 1
  Gl <- list(length = g.comp)
  for (i in 1:g.comp) {
    Gl[[i]] <- Reduce('+', mapply(FUN = "*", g[[i]], la.g[np.s.g[i]:np.e.g[i]]))
  }
  G <- spam::bdiag.spam(Matrix:::bdiag(Gl))
  return(G)
}

#' Helper function for ....
#'
#' @noRd
#' @keywords internal
construct.Ginv.matrix <- function(g,
                                  la) {
  la.g <- la[-1]
  g.comp <- length(g)
  g.tot.comp <- unlist(lapply(X = g, FUN = length))
  np.e.g <- cumsum(g.tot.comp)
  np.s.g <- np.e.g - g.tot.comp + 1
  Gl <- list(length = g.comp)
  for (i in 1:g.comp) {
    Gl[[i]] <- solve(Reduce('+', mapply(FUN = "*", g[[i]], la.g[np.s.g[i]:np.e.g[i]])))
  }
  Ginv <- spam::bdiag.spam(Matrix:::bdiag(Gl))
  return(Ginv)
}

#' Helper function for ....
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
