#' Compute tuning performance on training/validation data
#'
#' @param fit
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Y_list_train
#' @param indices_list_train
#'
#' @return
#' @export
compute_tuning_performance = function(fit, Y_list, X_list, indices_list, Y_list_train, indices_list_train) {

  n_gamma = fit$n_gamma
  n_lambda = fit$n_lambda

  p = nrow(fit[[1]][[1]]$Beta)
  q = ncol(fit[[1]][[1]]$Beta)

  if (ncol(X_list[[1]]) + 1 == nrow(fit$model_fits[[1]][[1]]$Beta)) {
    X_list = lapply(X_list, function(x) cbind(1, x))
  }

  n_iter = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  SSE = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  avg_R2 = matrix(NA, ncol = n_lambda, nrow = n_gamma)
  avg_correlation = matrix(NA, ncol = n_lambda, nrow = n_gamma)

  for (solution_path in fit$model_fits) {

    for (model in solution_path) {

      lambda = model$lambda_index
      gamma = model$gamma_index

      Beta = as.matrix(model$Beta)

      n_iter[gamma, lambda] = model$n_iter
      SSE[gamma, lambda] = compute_error(Y_list, X_list, indices_list, Beta)
      avg_R2[gamma, lambda] = compute_avg_R2(Y_list, X_list, indices_list, Y_list_train, indices_list_train, Beta)
      avg_correlation[gamma, lambda] = compute_avg_correlation(Y_list, X_list, indices_list, Beta)

    }

  }

  return(list(n_iter = n_iter,
              SSE = SSE,
              avg_R2 = avg_R2,
              avg_correlation = avg_correlation))

}

compute_gamma_weights = function(Y_list) {

  weights = sapply(Y_list, function(k) svd(k)$d[1])

  return(weights)

}

compute_candidate_gamma_sequence = function(n_gamma, min_ratio) {

  gamma = log_seq(1, min_ratio, n_gamma)

  return(gamma)

}

compute_candidate_lambda_grid = function(Y_list, X_list, q, indices_list, XtY_list, gamma_weights, gamma_sequence, n_lambda, min_ratio) {

  lambda = t(sapply(gamma_sequence, function(gamma) compute_candidate_lambda_sequence_fixed_gamma(Y_list, X_list, q, indices_list, XtY_list, n_lambda, min_ratio, gamma_weights, gamma)))

  return(lambda)

}

compute_candidate_lambda_sequence_fixed_gamma = function(Y_list, X_list, q, indices_list, XtY_list, n_lambda, min_ratio, gamma_weights, gamma) {

  result = matrix(0, nrow = ncol(X_list[[1]]), ncol = q)

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    L_k = nuclear_prox(Y_list[[k]], gamma_k)

    result[, indices_list[[k]]] = result[, indices_list[[k]]] + crossprod(X_list[[k]], L_k) - XtY_list[[k]]

  }

  max_lambda = max(abs(result[-1, ]))

  lambda = log_seq(max_lambda, max_lambda * min_ratio, n_lambda)

  return(lambda)

}

compute_candidate_lambda_sequence_glmnet = function(Y_list, X_list, q, indices_list, n_lambda, min_ratio) {

  result = matrix(0, nrow = ncol(X_list[[1]]), ncol = q)

  for (k in 1:length(Y_list)) {

    result[, indices_list[[k]]] = result[, indices_list[[k]]] + crossprod(X_list[[k]], Y_list[[k]])

  }

  max_lambda = max(abs(result))

  lambda = log_seq(max_lambda, max_lambda * min_ratio, n_lambda)

  return(lambda)

}

log_seq = function(from, to, length) {

  exp(seq(log(from), log(to), length.out = length))

}
