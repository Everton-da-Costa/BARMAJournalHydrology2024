#' @title Compute Log-Likelihood Terms for Beta AR Model (Auxiliary)
#' @description
#' Calculates the individual log-likelihood terms for a Beta Autoregressive
#' (BAR) model. This is an internal helper function.
#' These terms can be summed for the total log-likelihood.
#'
#' @param link A string specifying the link function (e.g., "logit", "probit",
#'             "cloglog", "loglog"). This defines the transformation
#'             \eqn{\eta = g(\mu)} linking the mean response \eqn{\mu \in (0,1)}
#'             to the linear predictor \eqn{\eta}. It's used to derive the link
#'             function \eqn{g(.)}, the inverse link (mean function)
#'             \eqn{\mu = g^{-1}(\eta)}, and its derivative \eqn{d\mu/d\eta}.
#' @param ar Numeric vector, AR orders to be included in the model.
#' @param alpha Scalar, the intercept parameter.
#' @param varphi Numeric vector, AR coefficients (\eqn{\varphi}) corresponding
#'               to lags in `ar`.
#' @param phi Scalar, the precision parameter (\eqn{\phi}).
#'
#' @details
#' This function takes the time series `y` (with values in the interval (0,1))
#' and a `link` specification (e.g., "logit", "probit", "cloglog") as input.
#' It internally generates `linkfun` (g(.)) and `linkinv` (the inverse of g(.),
#' also known as the mean function \eqn{\mu = g^{-1}(\eta)}) using the
#' provided `link` argument. The transformed time series `ynew = linkfun(y)`
#' and the `linkinv` function are then used to calculate the log-likelihood
#' terms.`error` are on the predictor scale.
#' This function is particularly useful for tests like LR2 in Costa et al.
#' (2024), which may require summing specific subsets of log-likelihood
#' terms from a null model.
#' @return A numeric vector of individual log-likelihood terms.
#' @keywords internal
loglik_terms_ar <- function(y, ar,link,
                            alpha = 0, varphi = 0, phi = 0) {

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
  n <- length(y)
  m <- max(p, na.rm = T)

  ynew <- linkfun(y)

  # --------------------------------------------------------------------- #
  eta <- rep(NA, n)
  mu <- rep(NA, n)

  for (t in (m + 1):n) {
    eta[t] <- alpha + varphi %*% ynew[t - ar]
    mu[t] <- linkinv(eta[t])
  }

  mu1 <- mu[(m + 1):n]
  y1 <- y[(m + 1):n]

  # --------------------------------------------------------------------- #
  # log-likelihood - function
  # --------------------------------------------------------------------- #

  ll_terms_ar <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

  return(ll_terms_ar)
}
