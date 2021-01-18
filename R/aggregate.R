aggregate_tuning_results = function(tuning_results, n_lambda, n_gamma) {

  p = nrow(tuning_results[[1]]$Beta)
  q = ncol(tuning_results[[1]]$Beta)

  aggregated_n_iter = matrix(NA, nrow = n_lambda, ncol = n_gamma)

  validation_error = matrix(NA, nrow = n_lambda, ncol = n_gamma)
  avg_validation_R2 = matrix(NA, nrow = n_lambda, ncol = n_gamma)
  weighted_avg_validation_R2 = matrix(NA, nrow = n_lambda, ncol = n_gamma)

  for (solution_path in tuning_results) {

    for (model in solution_path) {

      lambda = model$lambda_index
      gamma = model$gamma_index

      aggregated_n_iter[gamma, lambda] = model$n_iter

      if (!is.null(model$validation_error)) validation_error[gamma, lambda] = model$validation_error
      if (!is.null(model$avg_validation_R2)) avg_validation_R2[gamma, lambda] = model$avg_validation_R2
      if (!is.null(model$weighted_avg_validation_R2)) weighted_avg_validation_R2[gamma, lambda] = model$weighted_avg_validation_R2

    }

  }

  return(list(model_fits = tuning_results,
              validation_error = validation_error,
              avg_validation_R2 = avg_validation_R2,
              weighted_avg_validation_R2 = weighted_avg_validation_R2,
              n_iter = aggregated_n_iter))

}