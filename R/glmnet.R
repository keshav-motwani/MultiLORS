#' Fit glmnet model on concatenated data
#'
#' @param Y_list
#' @param X_list
#' @param D_list
#' @param Y_list_validation
#' @param X_list_validation
#' @param D_list_validation
#'
#' @return
#' @export
fit_glmnet = function(Y_list,
                      X_list,
                      indices_list,
                      Y_list_validation,
                      X_list_validation,
                      indices_list_validation,
                      n_lambda = 20,
                      lambda_min_ratio = 0.001) {

  p = ncol(X_list[[1]])
  q = max(unlist(indices_list))
  if (!is.null(X_list_validation)) {
    X_list_validation = lapply(X_list_validation, function(k) cbind(1, k))
  }

  lambda_sequence = compute_candidate_lambda_sequence_glmnet(Y_list, standardize_X(X_list), q, indices_list, n_lambda, lambda_min_ratio)

  Beta = array(dim = c(p + 1, q, n_lambda))

  for (index in 1:q) {

    data = subset_observed_data_univariate(Y_list, X_list, indices_list, index)

    fit = fit_glmnet_univariate(data$Y, data$X, lambda_sequence)
    Beta[, index, ] = fit

  }

  validation_error = numeric(n_lambda)
  avg_validation_R2 = numeric(n_lambda)
  weighted_avg_validation_R2 = numeric(n_lambda)

  models = list()

  for (lambda in 1:n_lambda) {

    fit = list(Beta = as(Beta[, , lambda], "dgCMatrix"))
    colnames(fit$Beta) = attr(indices_list, "responses")

    if (!is.null(Y_list_validation)) {

      error = compute_error(Y_list_validation, X_list_validation, indices_list_validation, Beta[, , lambda])
      avg_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, Beta[, , lambda])
      weighted_avg_R2 = compute_weighted_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, Beta[, , lambda])
      R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, Beta[, , lambda])

      validation_error[lambda] = error
      avg_validation_R2[lambda] = avg_R2
      weighted_avg_validation_R2[lambda] = weighted_avg_R2

      fit$validation_error = error
      fit$avg_validation_R2 = avg_R2
      fit$weighted_avg_validation_R2 = weighted_avg_R2
      fit$R2 = R2

    }

    models = c(models, list(fit))

  }

  return(list(model_fits = models,
              validation_error = validation_error,
              avg_validation_R2 = avg_validation_R2,
              weighted_avg_validation_R2 = weighted_avg_validation_R2,
              lambda_sequence = lambda_sequence))

}

#' @importFrom glmnet glmnet
fit_glmnet_univariate = function(Y, X, lambda_sequence) {

  fit = glmnet(X, Y, lambda = lambda_sequence/length(Y))

  return(as.matrix(coef(fit)))

}
