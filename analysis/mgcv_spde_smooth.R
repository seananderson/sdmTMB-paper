# -------------------------------------------------------------------------
# Code in this chunk is from:
#
# Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial
# Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020).
# https://doi.org/10.1007/s13253-019-00377-z
#
# Re-used here under a Creative Commons Attribution 4.0 International License
# http://creativecommons.org/licenses/by/4.0/
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  dim <- length(object$term)
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2)
    x[, 1] <- data[[object$term[1]]]
    x[, 2] <- data[[object$term[2]]]
  }
  if (is.null(object$xt)) {
    if (dim == 1) {
      t <- seq(min(x), max(x), len = object$bs.dim)
      mesh <- INLA::inla.mesh.1d(loc = t, degree = 2, boundary = "free")
    } else {
      stop("For 2D, mesh must be supplied as argument xt$mesh in s(...,xt = )")
    }
  } else {
    if (class(object$xt$mesh) != "inla.mesh") stop("xt must be NULL or an inla.mesh object")
    mesh <- object$xt$mesh
  }
  object$X <- as.matrix(INLA::inla.spde.make.A(mesh, x))
  inlamats <- INLA::inla.mesh.fem(mesh)
  object$S <- list()
  object$S[[1]] <- as.matrix(inlamats$c1)
  object$S[[2]] <- 2 * as.matrix(inlamats$g1)
  object$S[[3]] <- as.matrix(inlamats$g2)
  object$L <- matrix(c(2, 2, 2, 4, 2, 0), ncol = 2)
  object$rank <- rep(object$bs.dim, 3)
  object$null.space.dim <- 0
  object$mesh <- mesh
  object$df <- ncol(object$X)
  class(object) <- "spde.smooth"
  return(object)
}
Predict.matrix.spde.smooth <- function(object, data) {
  dim <- length(object$term)
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2)
    x[, 1] <- data[[object$term[1]]]
    x[, 2] <- data[[object$term[2]]]
  }
  Xp <- INLA::inla.spde.make.A(object$mesh, x)
  return(as.matrix(Xp))
}
# End of code from
# Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial
# Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020).
# https://doi.org/10.1007/s13253-019-00377-z
#
# Re-used here under a Creative Commons Attribution 4.0 International License
# http://creativecommons.org/licenses/by/4.0/
# -------------------------------------------------------------------------
