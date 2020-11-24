#' Aggregate LORS results computed with different lambda, gamma
#'
#' @param tuning_results
#' @param Y_list_validation
#' @param X_list_validation
#' @param D_list_validation
#'
#' @return
#' @export
aggregate_tuning_results = function(tuning_results) {

  n_lambda = max(sapply(tuning_results, `[[`, "lambda_index"))
  n_gamma = max(sapply(tuning_results, `[[`, "gamma_index"))
  p = nrow(tuning_results[[1]]$Beta)
  q = ncol(tuning_results[[1]]$Beta)

  aggregated_Beta = array(dim = c(n_lambda, n_gamma, p, q))
  aggregated_L_list = lapply(tuning_results[[1]]$L_list, function(k) array(dim = c(n_lambda, n_gamma, nrow(k), ncol(k))))
  aggregated_nuclear_norm_penalty = matrix(NA, nrow = n_lambda, ncol = n_gamma)
  aggregated_l1_penalty = matrix(NA, nrow = n_lambda, ncol = n_gamma)
  aggregated_n_iter = matrix(NA, nrow = n_lambda, ncol = n_gamma)
  validation_error = matrix(NA, nrow = n_lambda, ncol = n_gamma)

  for (model in tuning_results) {

    lambda = model$lambda_index
    gamma = model$gamma_index

    aggregated_Beta[lambda, gamma, 1:p, 1:q] = as.numeric(model$Beta)

    for (k in 1:length(model$L_list)) {
      aggregated_L_list[[k]][lambda, gamma, , ] = model$L_list[[k]]
    }

    aggregated_nuclear_norm_penalty[lambda, gamma] = nuclear_norm_penalty(model$L_list, 1, rep(1, length(model$L_list)))

    aggregated_l1_penalty[lambda, gamma] = l1_penalty(model$Beta[-1, ], 1)

    aggregated_n_iter[lambda, gamma] = model$n_iter

    validation_error[lambda, gamma] = model$validation_error

  }

  best_lambda_gamma = which(validation_error == min(validation_error, na.rm = TRUE), arr.ind = TRUE)
  best_Beta = aggregated_Beta[best_lambda_gamma[1], best_lambda_gamma[2], , ]

  print(best_lambda_gamma)

  return(list(Beta = aggregated_Beta,
              L_list = aggregated_L_list,
              nuclear_norm_penalty = aggregated_nuclear_norm_penalty,
              l1_penalty = aggregated_l1_penalty,
              n_iter = aggregated_n_iter,
              validation_error = validation_error,
              best_lambda = best_lambda_gamma[1],
              best_gamma = best_lambda_gamma[2],
              best_Beta = best_Beta))

}