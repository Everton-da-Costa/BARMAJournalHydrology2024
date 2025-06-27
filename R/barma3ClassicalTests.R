# =============================================================================
# Function for Classical Hypothesis Tests in \betaARMA Models
# =============================================================================

#' @title Classical Hypothesis Tests in \eqn{\beta}ARMA Models
#' @md
#'
#' @description
#' This R function, `barma3ClassicalTests`, performs a suite of classical
#' hypothesis tests -- Wald(W), Likelihood Ratio (LR), and Rao Score (RS) --
#' for Beta Autoregressive Moving Average (\eqn{\beta}ARMA) models. These
#' models are particularly useful for doubly-bounded time series, such as
#' those found in hydro-environmental studies (e.g., reservoir useful
#' volumes), which may exhibit seasonality and be subject to abnormal
#' regimes like prolonged droughts.
#'
#' The function is derived from an original BARMA implementation by Fabio M.
#' Bayer and has been significantly refactored and extended by Everton da Costa.
#'
#' @details
#' Key features and considerations based on recent research:
#' - Code structuring adhering to the tidyverse style guide for improved
#'   readability and maintainability.
#' - Implementation validated against the foundational paper:
#'   Rocha, A.V.; Cribari-Neto, F. (2009). Beta autoregressive moving
#'   average models. TEST, 18(3), 529--545. (Erratum: TEST, 26, 2017,
#'   451--459).
#' - This function implements the classical hypothesis tests presented in the
#'   study: Costa, E.; Cribari-Neto, F.; Scher, V.T. (2024). Test inferences 
#'   and link function selection in dynamic beta modeling of seasonal
#'   hydro-environmental time series with temporary abnormal regimes.
#'   Journal of Hydrology, 638, 131489.
#'   The study highlights potential inaccuracies in standard LR and Score tests
#'   when the number of individual conditional log-likelihoods differs
#'   between the null and non-null (unrestricted) models. This issue arises
#'   when autoregressive (p) or moving average (q) orders change, affecting
#'   `a = max(p,q)`. To address this, the function provides several test
#'   variants, including the recommended `LR_mfun_res`, `R_mfun_res`, and
#'   `R_exp_mfun_res`, which ensure comparability and lead to more accurate
#'   inferences.
#' - Optimization of the likelihood estimation process using the L-BFGS
#'   algorithm via the `lbfgs` package. (Ensure `lbfgs` is listed in
#'   Imports in DESCRIPTION if this is a package).
#' - Improved computational efficiency by minimizing redundant calculations.
#' - Comprehensive implementation of Wald, Likelihood Ratio, and Rao Score
#'   tests, including several variations, to assess parameter
#'   restrictions (e.g., significance of AR or MA components).
#'
#' @section Authorship and Version:
#' - Based on the original code by: Fabio M. Bayer (bayer@ufsm.br),
#'   Date: 2015-10-15
#' - Modified and improved by: Everton da Costa (everto.cost@gmail.com),
#'   Date: 2022-02-03
#' - VERSION: 1.00
#' - LAST UPDATE: 2023-07-27
#'
#' @param y A time series object of class `ts`. Values must be strictly within
#'        the interval (0, 1).
#' @param ar (Optional) A scalar or vector specifying the autoregressive (AR)
#'         orders for the unrestricted model.
#'         Example: `1` for AR(1), `c(1,3)` for AR with lags 1 and 3.
#'         Use `NA` if no AR component is desired.
#' @param ma (Optional) A scalar or vector specifying the moving average (MA)
#'         orders for the unrestricted model.
#'         Example: `1` for MA(1), `c(1,2)` for MA with lags 1 and 2.
#'         Use `NA` if no MA component is desired.
#' @param link A string specifying the link function for the conditional mean.
#'           Supported: `"logit"` (default), `"probit"`,`"cloglog"`,`"loglog"`.
#' @param rest_ar A scalar or vector indicating specific AR lags
#'              (from `ar` of unrestricted model) to restrict to zero
#'              under the null hypothesis.
#'              Example: If `ar = c(1,2)` and `rest_ar = 2`, tests
#'              H0: \eqn{\varphi_2 = 0}.
#'              If `rest_ar = c(1,2)`, tests H0: \eqn{\varphi_1=\varphi_2 = 0}.
#'              Defaults to `NA`.
#' @param rest_ma A scalar or vector indicating specific MA lags
#'              (from `ma` of unrestricted model) to restrict to zero
#'              under the null hypothesis. Similar to `rest_ar`.
#'              Defaults to `NA`.
#'              Note: Typically, either `rest_ar` OR `rest_ma` is specified.
#' @param arma_only A logical value (boolean).
#'              If `FALSE` (default): Performs full procedure: estimates
#'              unrestricted and restricted models, computes all tests.
#'              If `TRUE`: Primarily performs optimization for the specified
#'              BARMA(p,q) model (unrestricted_model output) and SKIPS tests.
#'              Other model/test outputs (e.g., `tests`,
#'              `ar_restricted_model`, `ma_restricted_model`) will be `NULL`
#'              or minimal.
#'
#' @return A list object with the following named elements:
#' \describe{
#'   \item{`unrestricted_model`}{Object (list) with detailed results from
#'                          fitting the full BARMA(p,q) model. Includes
#'                          estimates (`coeff`), SEs, log-likelihood
#'                          (`loglik`), convergence status (`conv`),
#'                          vcov matrix (`arma_vcov`).}
#'   \item{`ar_restricted_model`}{Object (list) with results from fitting a
#'                           BARMA model where AR components (specified by
#'                           `rest_ar` or all AR) are restricted to zero.
#'                           Estimated if AR restrictions are tested. May
#'                           be `NULL`.}
#'   \item{`ma_restricted_model`}{Object (list) with results from fitting a
#'                           BARMA model where MA components (specified by
#'                           `rest_ma` or all MA) are restricted to zero.
#'                           Estimated if MA restrictions are tested. May
#'                           be `NULL`.}
#'   \item{`tests`}{A list with Wald, Likelihood Ratio (LR), and Rao Score
#'             (RS) test results. Populated if `arma_only = FALSE` and
#'             models converge. Sub-elements refer to different
#'             formulations, some addressing issues from Costa et al. (2024)
#'             regarding the number of terms in log-likelihood
#'             (`a = max(p,q)`).
#'             \itemize{
#'               \item `W_res`: Wald test (stat, p-val). `W1` in Costa et
#'                        al. (2024), uses info matrix from unrestricted
#'                        model.
#'               \item `W2_res`: Alt. Wald test (stat, p-val). `W2` in
#'                         Costa et al. (2024), uses info matrix from
#'                         restricted model estimated with equal likelihood
#'                         terms (`mat_vcov_mfun`).
#'               \item `LR_naive_res`: "Naive" LR test (stat, p-val,
#'                               loglik_restr). `LR1` in Costa et al. (2024).
#'                               Uses `loglik_restr_naive` from restricted
#'                               model estimated with its own `a_N`. Can be
#'                               inaccurate or negative if `a_N < a_NN`.
#'               \item `LR_m0_res`: LR test (stat, p-val, loglik_restr
#'                            `loglik_restr_m0`). `LR2` from Costa
#'                            et al. (2024) if `loglik_restr_m0` is an
#'                            adjusted sum.
#'               \item `LR_mfun_res`: Recommended LR test (stat, p-val,
#'                              loglik_restr `loglik_restr_mfun`). `LR3` in
#'                              Costa et al. (2024). `loglik_restr_mfun`
#'                              from restricted model using `a_NN` (from
#'                              unrestricted model) for equal terms in
#'                              log-likelihoods. Non-negative.
#'               \item `R_naive_res`: "Naive" Rao Score (reduced form)
#'                              (stat, p-val). `Sr*` in Costa et al. (2024).
#'                              Uses `score_vec_naive`, `mat_vcov_naive`
#'                              from restricted model with its own `a_N`.
#'               \item `R_expanded_res`: "Naive" Rao Score (extended)
#'                                 (stat, p-val). `Se*` in Costa et al.
#'                                 (2024). Uses `score_vec_naive`,
#'                                 `mat_vcov_naive` from restricted model
#'                                 with its own `a_N`.
#'               \item `R_mfun_res`: Recommended Score (reduced) (stat,
#'                             p-val). `Sr` in Costa et al. (2024). Uses
#'                             `score_vec_mfun`, `mat_vcov_mfun` from
#'                             restricted model using `a_NN`.
#'               \item `R_exp_mfun_res`: Recommended Score (extended)
#'                                 (stat, p-val). `Se` in Costa et al.
#'                                 (2024). Uses `score_vec_mfun`,
#'                                 `mat_vcov_mfun` from restricted model
#'                                 using `a_NN`.
#'             }
#'             (Each result typically includes test statistic and p-value.
#'             DF for chi-squared tests is `num_rest`).
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'y' is a time series object (ts) with values in (0,1)
#' # Example: Test significance of AR(1) and MA(1) components
#' # in a BARMA(1,1) model.
#'
#' # To test AR components (H0: varphi_1 = 0)
#' # Note: This requires the actual barma3ClassicalTests function and its
#' # dependencies to be available and working.
#' # results_test_ar <- barma3ClassicalTests(y = y, ar = 1, ma = 1,
#' #                                        link = "logit", rest_ar = 1)
#' # print(results_test_ar$tests)
#'
#' # To test MA components (H0: theta_1 = 0)
#' # results_test_ma <- barma3ClassicalTests(y = y, ar = 1, ma = 1,
#' #                                        link = "logit", rest_ma = 1)
#' # print(results_test_ma$tests)
#'
#' # To estimate model without running specific tests (if arma_only = TRUE)
#' # model_only <- barma3ClassicalTests(y = y, ar = 1, ma = 1,
#' #                                   link = "logit", arma_only = TRUE)
#' # print(model_only$unrestricted_model)
#' }
#'
#' @references
#' Costa, E., Cribari-Neto, F., & Scher, V. T. (2024). Test inferences and
#' link function selection in dynamic beta modeling of seasonal
#' hydro-environmental time series with temporary abnormal regimes.
#' *Journal of Hydrology*, *638*, 131489.
#' \doi{10.1016/j.jhydrol.2024.131489}.
#'
#' Rocha, A.V. & Cribari-Neto, F. (2009). Beta autoregressive moving
#' average models. *TEST*, *18*(3), 529--545. (Erratum: *TEST*, *26*,
#' 2017, 451--459).
#'
#' @section BibTeX Citation for Costa et al. (2024):
#' For users who wish to cite the primary methodological paper by
#' Costa et al. (2024) that underpins aspects of this function, the
#' BibTeX entry is provided below for convenience:
#' ```
#' @article{Costa+Cribari+Scher_2024,
#'   title     = {Test inferences and link function selection in dynamic beta
#'                 modeling of seasonal hydro-environmental time series with
#'                 temporary abnormal regimes},
#'   author    = {Costa, E. and Cribari-Neto, F. and Scher, V. T.},
#'   journal   = {Journal of Hydrology},
#'   volume    = {638},
#'   pages     = {131489},
#'   year      = {2024},
#'   doi       = {10.1016/j.jhydrol.2024.131489}
#' }
#' ```
#' @export
#' @keywords \eqn{\beta}ARMA model; climate change; drought;
#' hydro-environmental data; link function; seasonality
#'
# =============================================================================


barma3ClassicalTests <- function(y, ar = NA, ma = NA, link = "logit",
                                 rest_ar = NA, rest_ma = NA, arma_only = FALSE) {
  # Link functions ============================================================

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

  maxit1 <- 1000

  # link functions - end ======================================================

  if (min(y) <= 0 || max(y) >= 1) stop("OUT OF RANGE (0,1)!")

  if (is.ts(y) == T) {
    freq <- frequency(y)
  } else {
    stop("data can be a time-series object")
  }

  if (any(is.na(ar)) == F) names_varphi <- c(paste("varphi", ar, sep = ""))
  if (any(is.na(ma)) == F) names_theta <- c(paste("theta", ma, sep = ""))

  p <- max(ar)
  q <- max(ma)
  n <- length(y)

  p1 <- length(ar)
  q1 <- length(ma)

  ynew <- linkfun(y)
  ystar <- log(y / (1 - y))

  # =========================================================================
  # Initializing the parameter values
  # =========================================================================
  # Start values for ARMA and BAR models
  start_value <- function(y) {
    # from:
    #       Beta Regression for Modelling Rates and Proportions
    #       Silvia Ferrari & Francisco Cribari-Neto
    #       p. 805

    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    x_start <- cbind(x_inter, P)
    y_start <- linkfun(y[(m + 1):n])

    fit_start <- lm.fit(x = x_start, y = y_start)

    mqo <- fit_start$coef

    alpha_start <- mqo[1]
    varphi_start <- mqo[-1]

    # --------------------------------------------------- #
    # precision
    # --------------------------------------------------- #
    k <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # ------------------------------------------------------------------------

    theta_start <- rep(0, q1)

    values_unrest <- c(alpha_start, varphi_start, theta_start, phi_start)
    values_rest <- c(alpha_start, varphi_start, phi_start)

    # ------------------------------------------------------------------------

    values_final <- list(
      values_unrest = values_unrest,
      values_rest = values_rest
    )

    return(values_final)
  }

  # out put ----------------------------------------------------------------- #
  unrestricted_model <- list()
  ar_restricted_model <- list()
  ma_restricted_model <- list()

  # z_info_K <- list()
  z_test <- list()

  # ========================================================================= #
  # ARMA model
  # ========================================================================= #

  if (any(is.na(ar) == F) && any(is.na(ma) == F)) {
    m <- max(p, q, na.rm = T)

    P <- matrix(nrow = n - m, ncol = p1)
    for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

    start_par <- start_value(y)
    # print(start_par$values_unrest)

    # ----------------------------------------------------------------------- #
    # log-likelihood
    # ----------------------------------------------------------------------- #
    loglik_arma <- function(z) {
      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta <- z[(p1 + 2):(p1 + q1 + 1)]
      phi <- z[p1 + q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
        mu[t] <- linkinv(eta[t])
        error[t] <- ynew[t] - eta[t]
      }

      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #

      ll_terms_arma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      final <- -1 * sum(ll_terms_arma)

      return(final)
    }

    # ----------------------------------------------------------------------- #
    # Vector Score
    # ----------------------------------------------------------------------- #
    score_arma <- function(z) {
      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      theta <- z[(p1 + 2):(p1 + q1 + 1)]
      phi <- z[p1 + q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
        mu[t] <- linkinv(eta[t])
        error[t] <- ynew[t] - eta[t] # preditor scale
        # error[t]   <- y[t] - mu[t]           # original scale
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]


      # --------------------------------------------------------------------- #
      mT <- diag(mu.eta(eta1))
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)

      R <- matrix(nrow = n - m, ncol = q1)
      for (t in 1:(n - m)) R[t, ] <- error[t + m - ma]

      # --------------------------------------------------------------------- #
      # Recorrencias
      # --------------------------------------------------------------------- #
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

      # --------------------------------------------------------------------- #
      ystar_mustar <- matrix(ystar - mustar, nrow = n - m, ncol = 1)

      U_alpha <- phi * s %*% mT %*% ystar_mustar
      U_varphi <- phi * t(rP) %*% mT %*% ystar_mustar
      U_theta <- phi * t(rR) %*% mT %*% ystar_mustar

      U_phi <- sum(mu1 * ystar_mustar + log(1 - y1)
        - digamma((1 - mu1) * phi) + digamma(phi))

      escore_vec <- c(U_alpha, U_varphi, U_theta, U_phi)

      # --------------------------------------------------------------------- #

      final <- -1 * escore_vec

      return(final)
    }

    # -------------------------------------------------------------------------

    names_par_arma <- c("alpha", names_varphi, names_theta, "phi")

    # ----------------------------------------------------------------------- #
    # optimization
    # ----------------------------------------------------------------------- #

    opt_arma <- lbfgs::lbfgs(
      call_eval = loglik_arma,
      call_grad = score_arma,
      vars = start_par$values_unrest,
      invisible = 1,
      # epsilon = 1e-6,
      # linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      max_iterations = maxit1
    )

    unrestricted_model$conv <- opt_arma$convergence

    # check convergence  --------------------------------------------------- #
    if (unrestricted_model$conv != 0) {
      warning("ARMA - FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      unrestricted_model$inv_inf_matrix <- 2
      # code 2 implies that don't have matrix to the tests
    } else if (unrestricted_model$conv == 0) {
      coeff_arma <- opt_arma$par[1:(p1 + q1 + 2)]
      names(coeff_arma) <- names_par_arma

      unrestricted_model$coeff <- coeff_arma
      unrestricted_model$loglik <- -1 * opt_arma$value

      # ------------------------------------------------------------------------
      # Information Fisher Matrix
      # ------------------------------------------------------------------------
      alpha <- coeff_arma[1]
      varphi <- coeff_arma[2:(p1 + 1)]
      theta <- coeff_arma[(p1 + 2):(p1 + q1 + 1)]
      phi <- coeff_arma[p1 + q1 + 2]

      # --------------------------------------------------------------------- #
      unrestricted_model$alpha <- alpha
      unrestricted_model$varphi <- varphi
      unrestricted_model$theta <- theta
      unrestricted_model$phi <- phi

      arma_K <- inf_matrix_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = unrestricted_model$alpha,
        varphi = unrestricted_model$varphi,
        theta = unrestricted_model$theta,
        phi = unrestricted_model$phi
      )

      unrestricted_model$arma_K <- arma_K

      # ========================================================================
      # Model presentation
      # ========================================================================
      solve_matrix_arma <- try(solve(arma_K, tol = 1e-20), silent = TRUE)

      # check inv matrix ---------------------------------------------------- #
      if (!(typeof(solve_matrix_arma) == "double")) {
        warning("ARMA, FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
        unrestricted_model$inv_inf_matrix <- 1
      } else {
        unrestricted_model$inv_inf_matrix <- 0

        # arma_vcov <- chol2inv(chol(arma_K))
        arma_vcov <- solve(arma_K, tol = 1e-20)

        unrestricted_model$arma_vcov <- arma_vcov

        stderror_arma <- sqrt(diag(arma_vcov))

        z_arma_stderror <- stderror_arma
        z_arma_zstat <- abs(unrestricted_model$coeff / stderror_arma)
        z_arma_pvalues <- 2 * (1 - pnorm(z_arma_zstat))

        unrestricted_model$counts <- as.numeric(opt_arma$counts[1])

        model_presentation_arma <- cbind(
          round(unrestricted_model$coeff, 4),
          round(z_arma_stderror, 4),
          round(z_arma_zstat, 4),
          round(z_arma_pvalues, 4)
        )

        colnames(model_presentation_arma) <- c(
          "Estimate",
          "Std. Error",
          "z value",
          "Pr(>|z|)"
        )

        unrestricted_model$model <- model_presentation_arma

        unrestricted_model$link <- link
      }
      # check inv matrix - end -------------------------------------------- #
    }
    # check convergence - end ----------------------------------------------- #
  }

  # ============================================================================
  # BAR Model
  # ============================================================================

  if (any(is.na(ar) == F) && any(is.na(rest_ma) == F) &&
    arma_only == FALSE) {
    m <- max(p, na.rm = T)
    num_rest <- length(rest_ma)

    P <- matrix(nrow = n - m, ncol = p1)
    for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]

    start_par <- start_value(y)
    # print(start_par$values_rest)

    # ----------------------------------------------------------------------- #
    # log-likelihood
    # ----------------------------------------------------------------------- #
    loglik_ar <- function(z) {
      m <- max(p, na.rm = T)

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]

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

      sum_ll <- -1 * sum(ll_terms_ar)


      return(sum_ll)
    }

    # ----------------------------------------------------------------------- #
    # score vector
    # ----------------------------------------------------------------------- #
    score_ar <- function(z) {
      m <- max(p, na.rm = T)

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar]
        mu[t] <- linkinv(eta[t])
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)
      mT <- diag(mu.eta(eta1))

      # --------------------------------------------------------------------- #
      Ualpha <- phi * sum(mu.eta(eta1) * (ystar - mustar))
      Uvarphi <- phi * t(P) %*% mT %*% (ystar - mustar)
      Uphi <- sum(mu1 * (ystar - mustar) + log(1 - y1)
        - digamma((1 - mu1) * phi) + digamma(phi))

      score_vec <- c(Ualpha, Uvarphi, Uphi)

      # --------------------------------------------------------------------- #

      final <- -1 * score_vec

      return(final)
    }

    names_par_ar <- c("alpha", names_varphi, "phi")

    # ----------------------------------------------------------------------- #
    # optim
    # ----------------------------------------------------------------------- #

    opt_ar <- lbfgs::lbfgs(
      call_eval = loglik_ar,
      call_grad = score_ar,
      vars = start_par$values_rest,
      invisible = 1,
      # epsilon = 1e-6,
      # linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      max_iterations = maxit1
    )

    # ======================================================================= #
    # BAR - Estimates for LRT: m inside function
    # ======================================================================= #
    # ----------------------------------------------------------------------- #
    # log-likelihood
    # ----------------------------------------------------------------------- #
    loglik_ar_mfun <- function(z) {
      m <- max(p, q, na.rm = T)

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]

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
      ll_terms_ar <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      sum_ll <- -1 * sum(ll_terms_ar)


      return(sum_ll)
    }

    # ----------------------------------------------------------------------- #
    # score vector
    # ----------------------------------------------------------------------- #
    score_ar_mfun <- function(z) {
      m <- max(p, q, na.rm = T)

      P_mfun <- matrix(NA, nrow = n - m, ncol = p1)
      for (i in 1:(n - m)) P_mfun[i, ] <- ynew[i + m - ar]

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      varphi <- z[2:(p1 + 1)]
      phi <- z[p1 + 2]

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + varphi %*% ynew[t - ar]
        mu[t] <- linkinv(eta[t])
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)
      mT <- diag(mu.eta(eta1))

      # --------------------------------------------------------------------- #
      Ualpha <- phi * sum(mu.eta(eta1) * (ystar - mustar))
      Uvarphi <- phi * t(P_mfun) %*% mT %*% (ystar - mustar)
      Uphi <- sum(mu1 * (ystar - mustar) + log(1 - y1)
        - digamma((1 - mu1) * phi) + digamma(phi))

      score_vec <- c(Ualpha, Uvarphi, Uphi)

      # --------------------------------------------------------------------- #

      final <- -1 * score_vec

      return(final)
    }

    # ----------------------------------------------------------------------- #
    # optim
    # ----------------------------------------------------------------------- #

    opt_ar_mfun <- lbfgs::lbfgs(
      call_eval = loglik_ar_mfun,
      call_grad = score_ar_mfun,
      vars = start_par$values_rest,
      invisible = 1,
      # epsilon = 1e-6,
      # linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      max_iterations = maxit1
    )

    # ======================================================================= #

    ar_restricted_model$conv <- opt_ar$convergence
    ar_restricted_model$conv_mfun <- opt_ar_mfun$convergence

    # check convergence ----------------------------------------------------- #
    if (ar_restricted_model$conv != 0 | ar_restricted_model$conv_mfun != 0) {
      warning("BAR - FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      ar_restricted_model$inv_inf_matrix <- 2
      inv_inf_matrix_rest <- ar_restricted_model$inv_inf_matrix
      # code 2 implies that don't have matrix to the tests
    } else if (ar_restricted_model$conv == 0 & ar_restricted_model$conv_mfun == 0) {
      coeff_ar <- opt_ar$par[1:(p1 + 2)]
      coeff_ar_mfun <- opt_ar_mfun$par[1:(p1 + 2)]

      names(coeff_ar) <- names_par_ar
      names(coeff_ar_mfun) <- names_par_ar

      ar_restricted_model$coeff <- coeff_ar
      ar_restricted_model$coeff_mfun <- coeff_ar_mfun

      ar_restricted_model$loglik <- -1 * opt_ar$value
      ar_restricted_model$loglik_mfun <- -1 * opt_ar_mfun$value

      alpha <- coeff_ar[1]
      varphi <- coeff_ar[2:(p1 + 1)]
      phi <- coeff_ar[p1 + 2]

      ar_alpha_mfun <- coeff_ar_mfun[1]
      ar_varphi_mfun <- coeff_ar_mfun[2:(p1 + 1)]
      ar_phi_mfun <- coeff_ar_mfun[p1 + 2]

      ar_alpha <- alpha
      ar_varphi <- varphi
      ar_phi <- phi

      # --------------------------------------------------------------------- #
      ar_restricted_model$alpha <- alpha
      ar_restricted_model$varphi <- varphi
      ar_restricted_model$phi <- phi

      # just check
      # score_ar_coeff <- score_ar(coeff_ar)
      # ar_restricted_model$score_ar_coeff <- score_ar_coeff

      # ------------------------------------------------------------------------
      # Fisher information matrix - BAR
      # ------------------------------------------------------------------------
      etahat <- rep(NA, n)
      muhat <- rep(NA, n)
      errorhat <- rep(0, n) # E(error)=0

      for (t in (m + 1):n) {
        etahat[t] <- alpha + varphi %*% ynew[t - ar]
        muhat[t] <- linkinv(etahat[t])
        # errorhat[t] <- y[t] - linkinv(etahat[t])    # original scale
        errorhat[t] <- ynew[t] - etahat[t] # predictor scale
      }
      etahat1 <- etahat[(m + 1):n]
      muhat1 <- muhat[(m + 1):n]
      y1 <- y[(m + 1):n]


      # --------------------------------------------------------------------- #
      vI <- as.vector(rep(1, n - m))

      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma((1 - muhat1) * phi)

      mT <- diag(mu.eta(etahat1))

      W <- diag(phi * (psi1 + psi2)) %*% mT^2
      vc <- phi * (psi1 * muhat1 - psi2 * (1 - muhat1))
      D <- diag(psi1 * (muhat1^2) + psi2 * (1 - muhat1)^2 - trigamma(phi))

      # --------------------------------------------------------------------- #
      t_P <- t(P)

      # ---------------------------------------- #
      K_a_a <- as.matrix(phi * sum(diag(W)))

      K_p_a <- phi * t_P %*% W %*% vI
      K_p_p <- phi * t_P %*% W %*% P

      K_a_phi <- vI %*% mT %*% vc
      K_p_phi <- t_P %*% mT %*% vc

      K_phi_a <- K_a_phi
      K_a_p <- t(K_p_a)
      K_phi_p <- t(K_p_phi)

      K_phi_phi <- sum(diag(D))

      ar_K <- rbind(
        cbind(K_a_a, K_a_p, K_a_phi),
        cbind(K_p_a, K_p_p, K_p_phi),
        cbind(K_phi_a, K_phi_p, K_phi_phi)
      )

      ar_restricted_model$ar_K <- ar_K

      # ======================================================================
      # Score vector ARMA - Auxiliary
      # ======================================================================
      ar_restricted_model$ar_score_vec_arma_naive <- score_vector_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ar_alpha,
        varphi = ar_varphi,
        theta = rep(0, num_rest),
        phi = ar_phi
      )

      ar_restricted_model$ar_score_vec_arma_mfun <- score_vector_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ar_alpha_mfun,
        varphi = ar_varphi_mfun,
        theta = rep(0, num_rest),
        phi = ar_phi_mfun
      )

      # ======================================================================
      # Fisher Information ARMA - Auxiliary
      # ======================================================================

      ar_restricted_model$ar_inf_matrix_arma_naive <- inf_matrix_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ar_alpha,
        varphi = ar_varphi,
        theta = rep(0, num_rest),
        phi = ar_phi
      )

      ar_restricted_model$ar_inf_matrix_arma_mfun <- inf_matrix_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ar_alpha_mfun,
        varphi = ar_varphi_mfun,
        theta = rep(0, num_rest),
        phi = ar_phi_mfun
      )

      # =======================================================================
      # Model presentation
      # =======================================================================
      solve_matrix_ar <- try(solve(ar_K, tol = 1e-20), silent = TRUE)

      # check matrix ----------------------------------------------------------
      if (!(typeof(solve_matrix_ar) == "double")) {
        warning("BAR, FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
        ar_restricted_model$inv_inf_matrix <- 1
        inv_inf_matrix_rest <- ar_restricted_model$inv_inf_matrix
      } else {
        ar_restricted_model$inv_inf_matrix <- 0
        inv_inf_matrix_rest <- ar_restricted_model$inv_inf_matrix

        # vcov_ar   <- chol2inv(chol(ar_K))
        vcov_ar <- solve(ar_K, tol = 1e-20)

        ar_restricted_model$vcov <- vcov_ar

        stderror_ar <- sqrt(diag(vcov_ar))

        ar_stderror <- stderror_ar
        ar_zstat <- abs(ar_restricted_model$coeff / stderror_ar)
        ar_pvalues <- 2 * (1 - pnorm(ar_zstat))

        ar_restricted_model$counts <- as.numeric(opt_ar$counts[1])

        model_presentation_ar <- cbind(
          round(ar_restricted_model$coeff, 4),
          round(ar_stderror, 4),
          round(ar_zstat, 4),
          round(ar_pvalues, 4)
        )

        colnames(model_presentation_ar) <- c(
          "Estimate",
          "Std. Error",
          "z value",
          "Pr(>|z|)"
        )

        ar_restricted_model$model <- model_presentation_ar

        ar_restricted_model$link <- link

        # Objects 1 for tests =============================================== #
        rest <- rest_ma

        loglik_restr_naive <- ar_restricted_model$loglik
        loglik_restr_mfun <- ar_restricted_model$loglik_mfun

        ll_terms_ar <- loglik_terms_ar(
          ar = ar,
          y = y,
          link = link,
          alpha = ar_alpha,
          varphi = ar_varphi,
          phi = ar_phi
        )

        loglik_restr_m0 <- sum(ll_terms_ar[num_rest:(n - 1)])

        # print(ll_terms_ar[num_rest:(n - 1)])

        score_vec_naive <- ar_restricted_model$ar_score_vec_arma_naive
        score_vec_mfun <- ar_restricted_model$ar_score_vec_arma_mfun

        matrix_inf_naive <- ar_restricted_model$ar_inf_matrix_arma_naive
        matrix_inf_mfun <- ar_restricted_model$ar_inf_matrix_arma_mfun
      }
      # check inv matrix - end ---------------------------------------------- #
    }
    # check convergence - end ----------------------------------------------- #
    conv_rest <- ar_restricted_model$conv
  }

  # ============================================================================
  # MA model
  # ============================================================================

  if (any(is.na(ma) == F) && any(is.na(rest_ar) == F) &&
    arma_only == FALSE) {
    m <- max(q, na.rm = T)
    num_rest <- length(rest_ar)

    # Initializing the parameter values MA - new --------------------------- #
    x_inter <- matrix(1, nrow = n - m, ncol = 1)
    y_start <- linkfun(y[(m + 1):n])

    fit_start <- lm.fit(x = x_inter, y = y_start)

    mqo <- fit_start$coef

    alpha_start <- abs(mqo[1])

    # --------------------------------------------------- #
    # precision
    # --------------------------------------------------- #
    k <- length(mqo)
    n1 <- n - m

    y_hat_fit_start <- fitted(fit_start)
    mean_fit_start <- linkinv(y_hat_fit_start)

    linkfun_deriv_aux <- mu.eta(eta = linkfun(mu = mean_fit_start))
    linkfun_deriv <- 1 / linkfun_deriv_aux

    er <- residuals(fit_start)
    sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)

    phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
    phi_start <- phi_start_aux / n1

    # ------------------------------------------------------------------------

    theta_start <- rep(0, q1)

    values_rest <- c(alpha_start, theta_start, phi_start)

    # ------------------------------------------------------------------------

    values_final <- list(values_rest = values_rest)

    # ---------------------------------------------------------------------- #
    # log-likelihood - MA
    # ---------------------------------------------------------------------- #
    loglik_ma <- function(z) {
      # --------------------------------------------------------------------- #
      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2]

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


      # log-likelihood - function ------------------------------------------- #
      ll_terms_ma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)


      sum_ll <- -1 * sum(ll_terms_ma)

      return(sum_ll)
    }

    # --------------------------------------------------------------------- #
    # score vector - MA
    # --------------------------------------------------------------------- #
    score_ma <- function(z) {
      # --------------------------------------------------------------------- #
      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2] # precision parameter

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)
      error <- rep(0, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + theta %*% error[t - ma]
        mu[t] <- linkinv(eta[t])
        error[t] <- ynew[t] - eta[t] # preditor scale
        # error[t] <- y[t] - linkinv(eta[t])      # original scale
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]


      # --------------------------------------------------------------------- #
      mT <- diag(mu.eta(eta1))
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)


      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m)) R[i, ] <- error[i + m - ma]

      # --------------------------------------------------------------------- #
      # Recorrencias
      # --------------------------------------------------------------------- #
      deta.dalpha <- rep(0, n)
      deta.dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i - ma]
        deta.dtheta[i, ] <- R[(i - m), ] - theta %*% deta.dtheta[i - ma, ]
      }

      s <- deta.dalpha[(m + 1):n]
      rR <- deta.dtheta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      Ualpha <- phi * s %*% mT %*% (ystar - mustar)
      Utheta <- phi * t(rR) %*% mT %*% (ystar - mustar)
      Uphi <- sum(mu1 * (ystar - mustar) + log(1 - y1)
        - digamma((1 - mu1) * phi) + digamma(phi))


      escore_vec <- -1 * c(Ualpha, Utheta, Uphi)

      # --------------------------------------------------------------------- #

      return(escore_vec)
    }

    names_par_ma <- c("alpha", names_theta, "phi")

    # ----------------------------------------------------------------------- #
    # optim
    # ----------------------------------------------------------------------- #

    opt_ma <- lbfgs::lbfgs(
      call_eval = loglik_ma,
      call_grad = score_ma,
      vars = start_par$values_rest,
      invisible = 1,
      # epsilon = 1e-6,
      # linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      max_iterations = maxit1
    )

    # ======================================================================= #
    # MA - Estimates for LRT: m inside function
    # ======================================================================= #
    # ---------------------------------------------------------------------- #
    # log-likelihood - MA
    # ---------------------------------------------------------------------- #
    loglik_ma_mfun <- function(z) {
      m <- max(p, q, na.rm = T)

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2]

      # --------------------------------------------------------------------- #
      error <- rep(0, n)
      eta <- rep(NA, n)
      mu <- rep(NA, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + theta %*% error[t - ma]
        mu[t] <- linkinv(eta[t])
        error[t] <- ynew[t] - eta[t]
      }

      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]


      # log-likelihood - function ------------------------------------------- #
      ll_terms_ma <- dbeta(y1, mu1 * phi, (1 - mu1) * phi, log = TRUE)

      sum_ll <- -1 * sum(ll_terms_ma)

      return(sum_ll)
    }

    # --------------------------------------------------------------------- #
    # score vector - MA
    # --------------------------------------------------------------------- #
    score_ma_mfun <- function(z) {
      m <- max(p, q, na.rm = T)

      # --------------------------------------------------------------------- #
      alpha <- z[1]
      theta <- z[2:(q1 + 1)]
      phi <- z[q1 + 2] # precision parameter

      # --------------------------------------------------------------------- #
      eta <- rep(NA, n)
      mu <- rep(NA, n)
      error <- rep(0, n)

      for (t in (m + 1):n) {
        eta[t] <- alpha + theta %*% error[t - ma]
        mu[t] <- linkinv(eta[t])
        error[t] <- ynew[t] - eta[t] # preditor scale
        # error[t] <- y[t] - linkinv(eta[t])      # original scale
      }

      eta1 <- eta[(m + 1):n]
      mu1 <- mu[(m + 1):n]
      y1 <- y[(m + 1):n]


      # --------------------------------------------------------------------- #
      mT <- diag(mu.eta(eta1))
      ystar <- log(y1 / (1 - y1))
      mustar <- digamma(mu1 * phi) - digamma((1 - mu1) * phi)


      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m)) R[i, ] <- error[i + m - ma]

      # --------------------------------------------------------------------- #
      # Recorrencias
      # --------------------------------------------------------------------- #
      deta.dalpha <- rep(0, n)
      deta.dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i - ma]
        deta.dtheta[i, ] <- R[(i - m), ] - theta %*% deta.dtheta[i - ma, ]
      }

      s <- deta.dalpha[(m + 1):n]
      rR <- deta.dtheta[(m + 1):n, ]

      # --------------------------------------------------------------------- #
      Ualpha <- phi * s %*% mT %*% (ystar - mustar)
      Utheta <- phi * t(rR) %*% mT %*% (ystar - mustar)
      Uphi <- sum(mu1 * (ystar - mustar) + log(1 - y1)
        - digamma((1 - mu1) * phi) + digamma(phi))


      escore_vec <- -1 * c(Ualpha, Utheta, Uphi)

      # --------------------------------------------------------------------- #

      return(escore_vec)
    }

    # ----------------------------------------------------------------------- #
    # optim
    # ----------------------------------------------------------------------- #

    opt_ma_mfun <- lbfgs::lbfgs(
      call_eval = loglik_ma_mfun,
      call_grad = score_ma_mfun,
      vars = start_par$values_rest,
      invisible = 1,
      # epsilon = 1e-6,
      # linesearch_algorithm = "LBFGS_LINESEARCH_BACKTRACKING_ARMIJO",
      max_iterations = maxit1
    )

    # =========================================================================


    ma_restricted_model$conv <- opt_ma$convergence
    ma_restricted_model$conv_mfun <- opt_ma_mfun$convergence

    # check convergence ----------------------------------------------------- #
    if (opt_ma$conv != 0 | ma_restricted_model$conv_mfun != 0) {
      warning("MA - FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      ma_restricted_model$inv_inf_matrix <- 2
      inv_inf_matrix_rest <- ma_restricted_model$inv_inf_matrix
      # code 2 implies that don't have matrix to the tests
    } else if (ma_restricted_model$conv == 0 & ma_restricted_model$conv_mfun == 0) {
      coeff_ma <- (opt_ma$par)[1:(q1 + 2)]
      coeff_ma_mfun <- (opt_ma_mfun$par)[1:(q1 + 2)]

      names(coeff_ma) <- names_par_ma
      names(coeff_ma_mfun) <- names_par_ma

      ma_restricted_model$coeff <- coeff_ma
      ma_restricted_model$coeff_mfun <- coeff_ma_mfun

      ma_restricted_model$loglik <- -1 * opt_ma$value
      ma_restricted_model$loglik_mfun <- -1 * opt_ma_mfun$value

      alpha <- coeff_ma[1]
      theta <- coeff_ma[2:(q1 + 1)]
      phi <- coeff_ma[q1 + 2]

      ma_alpha_mfun <- coeff_ma_mfun[1]
      ma_theta_mfun <- coeff_ma_mfun[2:(q1 + 1)]
      ma_phi_mfun <- coeff_ma_mfun[q1 + 2]

      ma_alpha <- alpha
      ma_theta <- theta
      ma_phi <- phi

      # --------------------------------------------------------------------- #
      ma_restricted_model$alpha <- alpha
      ma_restricted_model$theta <- theta
      ma_restricted_model$phi <- phi

      # ------------------------------------------------------------------------
      # Fisher Matrix Information - MA
      # -----------------------------------------------------------------------
      errorhat <- rep(0, n) # E(error)=0
      etahat <- rep(NA, n)
      muhat <- rep(NA, n)

      for (t in (m + 1):n) {
        etahat[t] <- alpha + theta %*% errorhat[t - ma]
        muhat[t] <- linkinv(etahat[t])
        # errorhat[t] <- y[t] - linkinv(etahat[t])  # original scale
        errorhat[t] <- ynew[t] - etahat[t] # predictor scale
      }
      etahat1 <- etahat[(m + 1):n]
      muhat1 <- muhat[(m + 1):n]
      y1 <- y[(m + 1):n]

      # --------------------------------------------------------------------- #
      R <- matrix(nrow = n - m, ncol = q1)
      for (i in 1:(n - m)) R[i, ] <- errorhat[i + m - ma]

      # --------------------------------------------------------------------- #
      # Recorrencias
      # --------------------------------------------------------------------- #
      deta.dalpha <- rep(0, n)
      deta.dtheta <- matrix(0, ncol = q1, nrow = n)

      for (i in (m + 1):n) {
        deta.dalpha[i] <- 1 - theta %*% deta.dalpha[i - ma]
        deta.dtheta[i, ] <- R[(i - m), ] - theta %*% deta.dtheta[i - ma, ]
      }

      s <- deta.dalpha[(m + 1):n]
      rR <- deta.dtheta[(m + 1):n, ]
      # --------------------------------------------------------------------- #

      psi1 <- trigamma(muhat1 * phi)
      psi2 <- trigamma((1 - muhat1) * phi)

      mT <- diag(mu.eta(eta = etahat1))

      W <- diag(c(phi * (psi1 + psi2))) %*% mT^2
      vc <- phi * (psi1 * muhat1 - psi2 * (1 - muhat1))
      D <- diag(c(psi1 * (muhat1^2) + psi2 * (1 - muhat1)^2 - trigamma(phi)))

      # --------------------------------------------------------------------- #
      t_s <- t(s)
      t_rR <- t(rR)

      # ---------------------------------------- #
      K_a_a <- phi * t_s %*% W %*% s

      K_t_a <- phi * t_rR %*% W %*% s
      K_t_t <- phi * t_rR %*% W %*% rR

      K_a_phi <- t_s %*% mT %*% vc
      K_t_phi <- t_rR %*% mT %*% vc

      K_a_t <- t(K_t_a)
      K_phi_a <- t(K_a_phi)
      K_phi_t <- t(K_t_phi)
      # ---------------------------------------- #

      K_phi_phi <- sum(diag(D))

      ma_K <- rbind(
        cbind(K_a_a, K_a_t, K_a_phi),
        cbind(K_t_a, K_t_t, K_t_phi),
        cbind(K_phi_a, K_phi_t, K_phi_phi)
      )

      ma_restricted_model$ma_K <- ma_K


      # ====================================================================== #
      # Score Vector ARMA - Auxiliary
      # ====================================================================== #
      ma_restricted_model$ma_score_vec_arma_naive <- score_vector_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ma_alpha,
        varphi = rep(0, num_rest),
        theta = ma_theta,
        phi = ma_phi
      )

      ma_restricted_model$ma_score_vec_arma_mfun <- score_vector_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ma_alpha_mfun,
        varphi = rep(0, num_rest),
        theta = ma_theta_mfun,
        phi = ma_phi_mfun
      )

      # ===================================================================== #
      # Fisher Information Matrix ARMA - Auxiliary
      # ===================================================================== #

      ma_restricted_model$ma_inf_matrix_arma_naive <- inf_matrix_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ma_alpha,
        varphi = rep(0, num_rest),
        theta = ma_theta,
        phi = ma_phi
      )

      ma_restricted_model$ma_inf_matrix_arma_mfun <- inf_matrix_arma(
        ar = ar,
        ma = ma,
        y = y,
        link = link,
        alpha = ma_alpha_mfun,
        varphi = rep(0, num_rest),
        theta = ma_theta_mfun,
        phi = ma_phi_mfun
      )

      # =======================================================================
      # Model presentation
      # =======================================================================
      solve_matrix_ma <- try(solve(ma_K, tol = 1e-20), silent = TRUE)

      # check inv matrix ---------------------------------------------------- #
      if (!(typeof(solve_matrix_ma) == "double")) {
        warning("MA, FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
        ma_restricted_model$inv_inf_matrix <- 1
        inv_inf_matrix_rest <- ma_restricted_model$inv_inf_matrix
      } else {
        ma_restricted_model$inv_inf_matrix <- 0
        inv_inf_matrix_rest <- ma_restricted_model$inv_inf_matrix

        # vcov_ma   <- chol2inv(chol(ma_K))
        vcov_ma <- solve(ma_K, tol = 1e-20)

        ma_restricted_model$vcov <- vcov_ma

        stderror_ma <- sqrt(diag(vcov_ma))

        ma_stderror <- stderror_ma
        ma_zstat <- abs(ma_restricted_model$coeff / stderror_ma)
        ma_pvalues <- 2 * (1 - pnorm(ma_zstat))

        ma_restricted_model$counts <- as.numeric(opt_ma$counts[1])

        model_presentation_ma <- cbind(
          round(ma_restricted_model$coeff, 4),
          round(ma_stderror, 4),
          round(ma_zstat, 4),
          round(ma_pvalues, 4)
        )

        colnames(model_presentation_ma) <- c(
          "Estimate",
          "Std. Error",
          "z value",
          "Pr(>|z|)"
        )

        ma_restricted_model$model <- model_presentation_ma

        ma_restricted_model$link <- link

        # Objects 1 for tests =============================================== #
        rest <- rest_ar

        loglik_restr_naive <- ma_restricted_model$loglik
        loglik_restr_mfun <- ma_restricted_model$loglik_mfun

        ll_terms_ma <- loglik_terms_ma(
          ma = ma,
          y = y,
          link = link,
          alpha = ma_alpha,
          theta = ma_theta,
          phi = ma_phi
        )

        loglik_restr_m0 <- sum(ll_terms_ma[num_rest:(n - 1)])

        score_vec_naive <- ma_restricted_model$ma_score_vec_arma_naive
        score_vec_mfun <- ma_restricted_model$ma_score_vec_arma_mfun

        matrix_inf_naive <- ma_restricted_model$ma_inf_matrix_arma_naive
        matrix_inf_mfun <- ma_restricted_model$ma_inf_matrix_arma_mfun
      }
      # check inv matrix - end ---------------------------------------------- #
    }
    # check convergence - end ----------------------------------------------- #
    conv_rest <- ma_restricted_model$conv
  }

  # ===========================================================================
  # Classical Hypothesis Tests: Wald, LR, and Rao Score
  # Methodological details based on Costa et al. (2024) J.Hydrol. 638, 131489
  # ===========================================================================

  if (arma_only == FALSE &&
    !is.null(unrestricted_model$arma_vcov) && inv_inf_matrix_rest == 0 &&
    unrestricted_model$conv == 0 && conv_rest == 0) {
    # Proceed with tests if: not 'arma_only', models converged, vcov ready.

    # --- Prepare components for hypothesis tests ---
    # Unrestricted estimates (for Wald)
    nu_hat <- unrestricted_model$coeff
    # Unrestricted log-likelihood (for LR)
    loglik_unrest <- unrestricted_model$loglik

    # Inv. info matrix (naive restricted): uses a_N. For Sr*, Se*.
    mat_vcov_naive <- solve(matrix_inf_naive, tol = 1e-20)
    # Inv. info matrix (mfun restricted): uses a_NN. For Sr, Se (recommended).
    mat_vcov_mfun <- solve(matrix_inf_mfun, tol = 1e-20)

    # --- Wald Test (W1 from Costa et al., 2024) ---
    # Uses variance-covariance matrix from the unrestricted model.
    W_expec_inf_inv <- unrestricted_model$arma_vcov
    W_solve_matrix <- solve(W_expec_inf_inv[rest, rest], tol = 1e-20)
    W_stat <- t(nu_hat[rest]) %*% W_solve_matrix %*% nu_hat[rest]
    W_pval <- pchisq(W_stat, df = num_rest, lower.tail = FALSE)
    W_res <- c(W_stat, W_pval)
    names(W_res) <- c("W_stat", "W_pval")
    z_test$W_res <- W_res

    # --- Wald Test (W2 from Costa et al., 2024) ---
    # Uses var-cov matrix from restricted model (estimated with a_NN via mfun).
    W2_expec_inf_inv <- mat_vcov_mfun
    W2_solve_matrix <- solve(W2_expec_inf_inv[rest, rest], tol = 1e-20)
    W2_stat <- t(nu_hat[rest]) %*% W2_solve_matrix %*% nu_hat[rest]
    W2_pval <- pchisq(W2_stat, df = num_rest, lower.tail = FALSE)
    W2_res <- c(W2_stat, W2_pval)
    names(W2_res) <- c("W2_stat", "W2_pval")
    z_test$W2_res <- W2_res

    # --- Likelihood Ratio Test (LR1 "naive" from Costa et al., 2024) ---
    # Uses loglik from restricted model estimated with its own a_N.
    # Can be inaccurate or negative if a_N (null model) < a_NN (non-null).
    LR_naive_stat <- 2 * (loglik_unrest - loglik_restr_naive)
    LR_naive_pval <- pchisq(LR_naive_stat, df = num_rest, lower.tail = FALSE)
    LR_naive_res <- c(LR_naive_stat, LR_naive_pval, loglik_restr_naive)
    names(LR_naive_res) <- c("LR_naive_stat", "LR_naive_pval", "LR_naive_loglik")
    z_test$LR_naive_res <- LR_naive_res

    # --- Likelihood Ratio Test (LR2 "adjusted" from Costa et al., 2024) ---
    # Restricted model (est. with a_N) loglik (loglik_restr_m0) is adjusted:
    # sums only the last n-a_NN individual log-likelihoods.
    LR_m0_stat <- 2 * (loglik_unrest - loglik_restr_m0)
    LR_m0_pval <- pchisq(LR_m0_stat, df = num_rest, lower.tail = FALSE)
    LR_m0_res <- c(LR_m0_stat, LR_m0_pval, loglik_restr_m0)
    names(LR_m0_res) <- c("LR_m0_stat", "LR_m0_pval", "LR_m0_loglik")
    z_test$LR_m0_res <- LR_m0_res

    # --- Likelihood Ratio Test (LR3 "recommended" from Costa et al., 2024) ---
    # Restricted model estimated using a_NN (via mfun) for loglik_restr_mfun.
    # Ensures same number of likelihood terms; non-negative statistic.
    LR_mfun_stat <- 2 * (loglik_unrest - loglik_restr_mfun)
    LR_mfun_pval <- pchisq(LR_mfun_stat, df = num_rest, lower.tail = FALSE)
    LR_mfun_res <- c(LR_mfun_stat, LR_mfun_pval, loglik_restr_mfun)
    names(LR_mfun_res) <- c("LR_mfun_stat", "LR_mfun_pval", "LR_mfun_loglik")
    z_test$LR_mfun_res <- LR_mfun_res

    # --- Rao Score Test (Sr* "naive reduced" from Costa et al., 2024) ---
    # Uses score (score_vec_naive) & var-cov (mat_vcov_naive) from restricted
    # model estimated with its own a_N. Can lead to size distortions.
    R_naive_stat <- (t(score_vec_naive[rest]) %*%
      mat_vcov_naive[rest, rest] %*%
      score_vec_naive[rest])
    R_naive_pval <- pchisq(R_naive_stat, df = num_rest, lower.tail = FALSE)
    R_naive_res <- c(R_naive_stat, R_naive_pval)
    names(R_naive_res) <- c("R_naive_stat", "R_naive_pval")
    z_test$R_naive_res <- R_naive_res

    # --- Rao Score Test (Se* "naive extended" from Costa et al., 2024) ---
    # Uses score (score_vec_naive) & var-cov (mat_vcov_naive) from restricted
    # model estimated with its own a_N. Original form by C.R. Rao.
    R_exp_stat <- t(score_vec_naive) %*% mat_vcov_naive %*% score_vec_naive
    R_exp_pval <- pchisq(R_exp_stat, df = num_rest, lower.tail = FALSE)
    R_exp_res <- c(R_exp_stat, R_exp_pval)
    names(R_exp_res) <- c("R_exp_stat", "R_exp_pval")
    z_test$R_expanded_res <- R_exp_res

    # --- Rao Score Test (Sr "recommended reduced" from Costa et al., 2024) ---
    # Uses score (score_vec_mfun) & var-cov (mat_vcov_mfun) from restricted
    # model estimated with a_NN (via mfun). More accurate.
    R_mfun_stat <- (t(score_vec_mfun[rest]) %*%
      mat_vcov_mfun[rest, rest] %*%
      score_vec_mfun[rest])
    R_mfun_pval <- pchisq(R_mfun_stat, df = num_rest, lower.tail = FALSE)
    R_mfun_res <- c(R_mfun_stat, R_mfun_pval)
    names(R_mfun_res) <- c("R_mfun_stat", "R_mfun_pval")
    z_test$R_mfun_res <- R_mfun_res

    # --- Rao Score Test (Se "recommended extended" from Costa et al., 2024) ---
    # Uses score (score_vec_mfun) & var-cov (mat_vcov_mfun) from restricted
    # model estimated with a_NN (via mfun). More accurate.
    R_exp_mfun_stat <- t(score_vec_mfun) %*% mat_vcov_mfun %*% score_vec_mfun
    R_exp_mfun_pval <- pchisq(R_exp_mfun_stat, df = num_rest, lower.tail = FALSE)
    R_exp_mfun_res <- c(R_exp_mfun_stat, R_exp_mfun_pval)
    names(R_exp_mfun_res) <- c("R_exp_mfun_stat", "R_exp_mfun_pval")
    z_test$R_exp_mfun_res <- R_exp_mfun_res
  }

  final <- list(
    unrestricted_model = unrestricted_model,
    ar_restricted_model = ar_restricted_model,
    ma_restricted_model = ma_restricted_model,
    classical_tests = z_test
  )

  return(final)
}
