#' Fit LORS model
#'
#' @param Y_list
#' @param X_list
#' @param D_list
#' @param Y_list_validation
#' @param X_list_validation
#' @param D_list_validation
#' @param standardize
#' @param n_lambda
#' @param n_gamma
#' @param lambda_min_ratio
#' @param gamma_min_ratio
#' @param n_iter
#' @param tolerance
#' @param line_search
#' @param verbose
#'
#' @return
#' @export
MultiLORS = function(Y_list,
                     X_list,
                     D_list,
                     Y_list_validation = NULL,
                     X_list_validation = NULL,
                     D_list_validation = NULL,
                     standardize = TRUE,
                     n_lambda = 20,
                     n_gamma = 20,
                     lambda_min_ratio = 0.001,
                     gamma_min_ratio = 0.001,
                     n_iter = 1000,
                     tolerance = 1e-6,
                     early_stopping = TRUE,
                     line_search = TRUE,
                     verbose = FALSE) {

  if (standardize) X_list = standardize_X(X_list)
  X_mean = attributes(X_list)$mean
  X_sd = attributes(X_list)$sd
  X_list = lapply(X_list, function(k) cbind(1, k))
  XtX_list = lapply(X_list, crossprod)
  Y_dot_list = mapply(zero_pad_matrix, Y_list, D_list)
  XtY_dot_list = mapply(crossprod, x = X_list, y = Y_dot_list, SIMPLIFY = FALSE)
  indices_list = lapply(D_list, function(k) which(diag(k) == 1))
  if (!is.null(X_list_validation)) {
    X_list_validation = lapply(X_list_validation, function(k) cbind(1, k))
    indices_list_validation = lapply(D_list_validation, function(k) which(diag(k) == 1))
  }

  gamma_weights = compute_gamma_weights(Y_list)
  gamma_sequence = compute_candidate_gamma_sequence(n_gamma, gamma_min_ratio)
  lambda_grid = compute_candidate_lambda_grid(Y_list, X_list, D_list, XtY_dot_list, gamma_weights, gamma_sequence, n_lambda, lambda_min_ratio)

  s_Beta = compute_s_Beta(XtX_list, D_list)

  result = list()

  for (gamma in 1:n_gamma) {

    Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = ncol(D_list[[1]]))
    L_list_old = lapply(Y_list, function(k) matrix(0, nrow = nrow(k), ncol = ncol(k)))

    for (lambda in 1:n_lambda) {

      print(paste0("gamma: ", gamma, "; lambda: ", lambda))

      print(system.time({
        model = fit(
          Y_list = Y_list,
          X_list = X_list,
          D_list = D_list,
          indices_list = indices_list,
          XtX_list = XtX_list,
          XtY_dot_list = XtY_dot_list,
          X_mean = X_mean,
          X_sd = X_sd,
          lambda = lambda_grid[gamma, lambda],
          gamma = gamma_sequence[gamma],
          gamma_weights = gamma_weights,
          Beta_old = Beta_old,
          L_list_old = L_list_old,
          s_Beta = s_Beta,
          n_iter = n_iter,
          tolerance = tolerance,
          line_search = line_search,
          verbose = verbose
        )
      }))

      model$lambda_index = lambda
      model$gamma_index = gamma

      if (!is.null(Y_list_validation)) {
        model$validation_error = evaluate_g(Y_list_validation, X_list_validation, rep(0, length(Y_list_validation)), indices_list_validation, model$Beta)
        print(paste0("Validation Error: ", model$validation_error))
      }

      result = c(result, list(model))

      if (early_stopping && !is.null(Y_list_validation)) {
        if (lambda > 5 && lambda > n_lambda/4 && model$validation_error > 1.05 * result[[length(result) - 1]][["validation_error"]]) {
          break
        }
      }

      Beta_old = model$Beta
      L_list_old = model$L_list

    }
  }

  result = aggregate_tuning_results(result)

  return(result)

}
