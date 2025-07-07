#' @title Compute Score Vector for BARMA Model (Auxiliary)
#' @description
#' Calculates the score vector (gradient of the log-likelihood) for a BARMA
#' model. This is an internal helper function.
#' The score vector is used in ML estimation and for Score tests.
#'
#' @param link A string specifying the link function (e.g., "logit", "probit",
#'             "cloglog", "loglog"). This defines the transformation
#'             \eqn{\eta = g(\mu)} linking the mean response \eqn{\mu \in (0,1)}
#'             to the linear predictor \eqn{\eta}. It's used to derive the link
#'             function \eqn{g(.)}, the inverse link (mean function)
#'             \eqn{\mu = g^{-1}(\eta)}, and its derivative \eqn{d\mu/d\eta}.
#' @param alpha Scalar, the intercept parameter.
#' @param varphi Numeric vector, AR coefficients (\eqn{\varphi}).
#' @param theta Numeric vector, MA coefficients (\eqn{\theta}).
#' @param phi Scalar, the precision parameter (\eqn{\phi}).
#'
#' @details
#' This function uses the time series `y` (values in (0,1)) and a `link`
#' string (e.g., "logit") to internally derive `linkfun` (g(.)), `linkinv`
#' (\eqn{\mu = g^{-1}(\eta)}), and `mu.eta` (\eqn{d\mu/d\eta}). The transformed
#' series `ynew = linkfun(y)` is then used.
#' It computes components of the Score Vector for parameters
#' \eqn{\alpha, \varphi, \theta, \phi}. `error` are on the
#' predictor scale. `ystar` is `logit(y1)` and `mustar` involves digamma
#' functions.
#'
#' @return A numeric vector representing the score for parameters
#'         (\eqn{\alpha}, \eqn{\varphi}, \eqn{\theta}, \eqn{\phi}).
#' @keywords internal
score_vector_arma <- function(y, ar, ma, link,
                              alpha = 0, varphi = 0, theta = 0, phi = 0) {

  # Link functions
  # ----------------------------------------------------------------------- #
  # Check if the link argument is a character string or an expression
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      linktemp <- eval(link)
    }
  }

  # Set up the link function based on the provided link argument
  if (linktemp == "logit") {

    # Use make.link function for logit link (efficient C implementation)
    stats <- make.link("logit")

  } else if (linktemp == "probit") {

    # Use make.link function for probit link
    stats <- make.link("probit")

  } else if (linktemp == "cloglog") {

    # Use make.link function for cloglog link
    stats <- make.link("cloglog")

  } else if (linktemp == "loglog") {

    # Manually define linkfun, linkinv, and mu.eta for loglog link
    stats <- list()

    stats$linkfun <- function(mu) -log(-log(mu))
    stats$linkinv <- function(eta) exp(-exp(-eta))
    stats$mu.eta <- function(eta) exp(-exp(-eta)) * exp(-eta)

  } else {
    # If the provided link is not supported, raise an error
    stop(
      paste0(
        "The '", linktemp, "' link function is not available.", "\n",
        "Available options are: ", "logit, loglog, cloglog, and probit."
      )
    )
  }

  # Create a list with the link function details
  link1 <- structure(list(
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta = stats$mu.eta
  ))

  # Assign link function components to separate variables
  linkfun <- link1$linkfun
  linkinv <- link1$linkinv
  mu.eta <- link1$mu.eta

  # --------------------------------------------------------------------- #
  p <- max(ar)
  p1 <- length(ar)

  q <- max(ma)
  q1 <- length(ma)

  n <- length(y)
  m <- max(p, q, na.rm = T)

  ynew <- linkfun(y)

  # ----------------------------------------------------------------------- #
  P <- matrix(nrow = n - m, ncol = p1)
  for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]


  # ----------------------------------------------------------------------- #
  error <- rep(0, n)
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  for (t in (m + 1):n) {
    eta[t] <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
    mu[t] <- linkinv(eta[t])
    error[t] <- ynew[t] - eta[t]
  }

  eta1 <- eta[(m + 1):n]
  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # ------------------------------------------------------------------------ #

  mT <- diag(mu.eta(eta1))
  ystar <- log(y1 / (1 - y1))
  mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

  R <- matrix(nrow = n - m, ncol = q1)
  for (t in 1:(n - m)) R[t, ] <- error[t + m - ma]

  # Recorrencias --------------------------------------------------------- #

  deta.dalpha <- rep(0, n)
  deta.dvarphi <- matrix(0, nrow = n, ncol = p1)
  deta.dtheta <- matrix(0, nrow = n, ncol = q1)

  for (t in (m + 1):n) {
    deta.dalpha[t] <- 1 - theta %*% deta.dalpha[t - ma]
    deta.dvarphi[t, ] <- P[t - m, ] - theta %*% deta.dvarphi[t - ma, ]
    deta.dtheta[t, ] <- R[t - m, ] - theta %*% deta.dtheta[t - ma, ]
  }

  s <- deta.dalpha[(m + 1):n]
  rP <- deta.dvarphi[(m + 1):n, ]
  rR <- deta.dtheta[(m + 1):n, ]

  # ---------------------------------------------------------------------- #
  U_alpha <- phi * s %*% mT %*% (ystar - mustar)
  U_varphi <- phi * t(rP) %*% mT %*% (ystar - mustar)
  U_theta <- phi * t(rR) %*% mT %*% (ystar - mustar)

  U_phi <- sum(mu1 * (ystar - mustar) + log(1 - y1)
               - digamma((1 - mu1) * phi) + digamma(phi))

  # ---------------------------------------------------------------------- #

  escore_vec <- c(U_alpha, U_varphi, U_theta, U_phi)

  return(escore_vec)
}
