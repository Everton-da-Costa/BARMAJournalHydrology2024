#' @title Compute Log-Likelihood Terms for Beta MA Model (Auxiliary)
#' @description
#' Calculates the individual log-likelihood terms for a Beta Moving Average
#' (BMA) model. This is an internal helper function.
#'
#' @param ma Numeric vector, MA orders to be included in the model.
#' @param alpha Scalar, the intercept parameter.
#' @param theta Numeric vector, MA coefficients (\eqn{\theta}) corresponding
#'              to lags in `ma`.
#' @param phi Scalar, the precision parameter (\eqn{\phi}).
#'
#' @details
#' #' This function takes the time series `y` (with values in the interval
#' (0,1)) and a `link` specification (e.g., "logit", "probit", "cloglog") as
#' input. It internally generates `linkfun` (g(.)) and `linkinv` (the inverse
#' of g(.), also known as the mean function \eqn{\mu = g^{-1}(\eta)}) using the
#' provided `link` argument. The transformed time series `ynew = linkfun(y)`
#' and the `linkinv` function are then used to calculate the log-likelihood
#' terms. `error` are on the predictor scale. This function is particularly
#' useful for tests like LR2 in Costa et al. (2024), which may require summing
#'  specific subsets of log-likelihood terms from a null model.
#'
#' @return A numeric vector of individual log-likelihood terms.
#' @keywords internal
loglik_terms_ma <- function(y, ma, link,
                            alpha = 0, theta = 0, phi = 0) {

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

  q <- max(ma)
  q1 <- length(ma)

  n <- length(y)
  m <- max(q, na.rm = T)

  ynew <- linkfun(y)
  # --------------------------------------------------------------------- #
  error <- rep(0, n)
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  for (t in (m + 1):n) {
    eta[t] <- alpha + theta %*% error[t - ma]
    mu[t] <- linkinv(eta[t])
    error[t] <- ynew[t] - eta[t] # preditor scale
    # error[t] <- y[t] - linkinv(eta[t])        # original scale
  }

  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # --------------------------------------------------------------------- #
  ll_terms_ma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  # sum_ll <- -1 * sum(ll_terms_ma)

  return(ll_terms_ma)
}
