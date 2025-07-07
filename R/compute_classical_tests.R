#' @title Compute Classical Hypothesis Tests
#' @description Internal function to compute Wald, LR, and Rao Score tests.
#'   Methodological details based on Costa et al. (2024) J.Hydrol. 638, 131489.
#' @param unrestricted_model The fitted unrestricted model object.
#' @param loglik_restr_naive Log-likelihood of the restricted model (using a_N).
#' @param loglik_restr_m0 Adjusted log-likelihood of the restricted model.
#' @param loglik_restr_mfun Log-likelihood of the restricted model (using a_NN).
#' @param score_vec_naive Score vector from the restricted model (using a_N).
#' @param score_vec_mfun Score vector from the restricted model (using a_NN).
#' @param matrix_inf_naive Information matrix from the restricted model
#'   (using a_N).
#' @param matrix_inf_mfun Information matrix from the restricted model
#'   (using a_NN).
#' @param rest A vector of indices for the restricted parameters.
#' @param num_rest The number of restricted parameters (degrees of freedom).
#' @return A list (`z_test`) containing the results of all computed tests.
#' @keywords internal
compute_classical_tests <- function(unrestricted_model,
                                    loglik_restr_naive,
                                    loglik_restr_m0,
                                    loglik_restr_mfun,
                                    score_vec_naive,
                                    score_vec_mfun,
                                    matrix_inf_naive,
                                    matrix_inf_mfun,
                                    rest,
                                    num_rest) {
  z_test <- list()

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

  return(z_test)
}
