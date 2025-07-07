#' @title Compute Fisher Information Matrix for BARMA Model (Auxiliary)
#' @description
#' Calculates the expected Fisher Information Matrix for a Beta Autoregressive
#' Moving Average (BARMA) model. This is an internal helper function.
#' The matrix is derived based on the methodology in Rocha & Cribari-Neto
#' (2009, 2017) and is essential for statistical inference, e.g., for
#' calculating standard errors or for constructing Wald and Score tests
#' as discussed in Costa et al. (2024).
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
#' It computes components of the Fisher Information Matrix for parameters
#' \eqn{\alpha, \varphi, \theta, \phi}. `errorhat` are on the
#' predictor scale.
#'
#' @return A matrix representing the Fisher Information Matrix for the
#'         parameters (\eqn{\alpha}, \eqn{\varphi}, \eqn{\theta}, \eqn{\phi}).
#'
#' @references
#' Rocha, A.V., Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. TEST, 18(3), 529--545.
#' Costa, E., Cribari-Neto, F., Scher, V.T. (2024). Test inferences and
#' link function selection in dynamic beta modeling... Journal of Hydrology,
#' 638, 131489.
#' @keywords internal
inf_matrix_arma <- function(y, ar, ma, link,
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
  P <- matrix(NA, nrow = n - m, ncol = p1)
  for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]


  # ----------------------------------------------------------------------- #
  errorhat <- rep(0, n)
  muhat <- rep(0, n)
  etahat <- rep(NA, n)

  for (t in (m + 1):n) {
    etahat[t] <- alpha + varphi %*% ynew[t - ar] + theta %*% errorhat[t - ma]
    muhat[t] <- linkinv(etahat[t])
    errorhat[t] <- ynew[t] - etahat[t]
  }

  etahat1 <- etahat[(m + 1):n]
  muhat1 <- muhat[(m + 1):n]
  y1 <- y[(m + 1):n]


  # ----------------------------------------------------------------------- #
  R <- matrix(nrow = n - m, ncol = q1)
  for (i in 1:(n - m)) R[i, ] <- errorhat[i + m - ma]


  # Recorrencias ---------------------------------------------------------- #
  deta.dalpha <- rep(0, n)
  deta.dvarphi <- matrix(0, ncol = p1, nrow = n)
  deta.dtheta <- matrix(0, ncol = q1, nrow = n)

  for (i in (m + 1):n) {
    deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i - ma]
    deta.dvarphi[i, ] <- P[(i - m), ] - theta %*% deta.dvarphi[i - ma, ]
    deta.dtheta[i, ] <- R[(i - m), ] - theta %*% deta.dtheta[i - ma, ]
  }


  s <- deta.dalpha[(m + 1):n]
  rP <- deta.dvarphi[(m + 1):n, ]
  rR <- deta.dtheta[(m + 1):n, ]

  # Recorrencias - end ---------------------------------------------------- #

  psi1 <- trigamma(muhat1 * phi)
  psi2 <- trigamma((1 - muhat1) * phi)

  mT <- diag(mu.eta(eta = etahat1))

  W <- diag(phi * (psi1 + psi2)) %*% mT^2
  vc <- phi * (psi1 * muhat1 - psi2 * (1 - muhat1))
  D <- diag(psi1 * muhat1^2 + psi2 * (1 - muhat1)^2 - trigamma(phi))

  # ----------------------------------------------------------------------- #
  t_s <- t(s)
  t_rP <- t(rP)
  t_rR <- t(rR)

  # ---------------------------------------- #
  K_a_a <- phi * t_s %*% W %*% s
  K_p_a <- phi * t_rP %*% W %*% s
  K_t_a <- phi * t_rR %*% W %*% s

  K_p_p <- phi * t_rP %*% W %*% rP
  K_p_t <- phi * t_rP %*% W %*% rR
  K_t_t <- phi * t_rR %*% W %*% rR

  K_a_phi <- t_s %*% mT %*% vc
  K_p_phi <- t_rP %*% mT %*% vc
  K_t_phi <- t_rR %*% mT %*% vc

  K_phi_phi <- sum(diag(D))
  # ---------------------------------------- #
  K_a_p <- t(K_p_a)
  K_t_p <- t(K_p_t)
  K_a_t <- t(K_t_a)

  K_phi_a <- K_a_phi
  K_phi_p <- t(K_p_phi)
  K_phi_t <- t(K_t_phi)

  # ----------------------------------------------------------------------- #

  arma_K <- rbind(
    cbind(K_a_a, K_a_p, K_a_t, K_a_phi),
    cbind(K_p_a, K_p_p, K_p_t, K_p_phi),
    cbind(K_t_a, K_t_p, K_t_t, K_t_phi),
    cbind(K_phi_a, K_phi_p, K_phi_t, K_phi_phi)
  )
}
