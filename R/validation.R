#' Compute performance on validation data
#'
#' @param fit
#' @param Y_list_validation
#' @param X_list_validation
#' @param indices_list_validation
#' @param Y_list_train
#' @param indices_list_train
#'
#' @return
#' @export
compute_validation_performance = function(fit, Y_list_validation, X_list_validation, indices_list_validation, Y_list_train, indices_list_train) {

  n_gamma = fit$n_gamma
  n_lambda = fit$n_lambda

  p = nrow(fit[[1]][[1]]$Beta)
  q = ncol(fit[[1]][[1]]$Beta)

  if (ncol(X_list_validation[[1]]) + 1 == nrow(fit$model_fits[[1]][[1]]$Beta)) {
    X_list_validation = lapply(X_list_validation, function(x) cbind(1, x))
  }

  n_iter = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  validation_error = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  avg_validation_R2 = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  weighted_avg_validation_R2 = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  avg_validation_correlation = matrix(NA, ncol = n_lambda, nrow = n_gamma)

  for (solution_path in fit$model_fits) {

    for (model in solution_path) {

      lambda = model$lambda_index
      gamma = model$gamma_index

      Beta = as.matrix(model$Beta)

      n_iter[gamma, lambda] = model$n_iter
      validation_error[gamma, lambda] = compute_error(Y_list_validation, X_list_validation, indices_list_validation, Beta)
      avg_validation_R2[gamma, lambda] = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list_train, indices_list_train, Beta)
      weighted_avg_validation_R2[gamma, lambda] = compute_weighted_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list_train, indices_list_train, Beta)
      avg_validation_correlation[gamma, lambda] = compute_avg_correlation(Y_list_validation, X_list_validation, indices_list_validation, Beta)

    }

  }

  return(list(n_iter = n_iter,
              validation_error = validation_error,
              avg_validation_R2 = avg_validation_R2,
              weighted_avg_validation_R2 = weighted_avg_validation_R2,
              avg_validation_correlation = avg_validation_correlation))

}
