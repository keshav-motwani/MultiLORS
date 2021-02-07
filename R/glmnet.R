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

  lambda_sequence = compute_candidate_lambda_sequence_glmnet(Y_list, standardize_X(X_list), q, indices_list, n_lambda, lambda_min_ratio)

  Beta = array(dim = c(p + 1, q, n_lambda))

  for (index in 1:q) {

    data = subset_observed_data_univariate(Y_list, X_list, indices_list, index)

    fit = fit_glmnet_univariate(data$Y, data$X, lambda_sequence)
    Beta[, index, ] = fit

  }

  X_list = lapply(X_list, function(k) cbind(1, k))
  if (!is.null(X_list_validation)) {
    X_list_validation = lapply(X_list_validation, function(k) cbind(1, k))
  }

  train_SSE = numeric(n_lambda)
  avg_train_R2 = numeric(n_lambda)
  avg_train_correlation = numeric(n_lambda)

  validation_SSE = numeric(n_lambda)
  avg_validation_R2 = numeric(n_lambda)
  avg_validation_correlation = numeric(n_lambda)

  models = list()

  for (lambda in 1:n_lambda) {

    fit = list(Beta = as(Beta[, , lambda], "dgCMatrix"), performance = list(train = list(), validation = list()))
    colnames(fit$Beta) = attr(indices_list, "responses")
    if (!is.null(colnames(X_list[[1]]))) rownames(fit$Beta) = c("intercept", colnames(X_list[[1]]))

    R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, fit$Beta)
    correlation = compute_correlation(Y_list, X_list, indices_list, fit$Beta)

    fit$performance$train$R2 = R2
    fit$performance$train$correlation = correlation

    error = compute_error(Y_list, X_list, indices_list, as.matrix(fit$Beta))
    avg_R2 = compute_avg_R2(Y_list, X_list, indices_list, Y_list, indices_list, fit$Beta)
    avg_correlation = compute_avg_correlation(Y_list, X_list, indices_list, fit$Beta)

    train_SSE[lambda] = error
    avg_train_R2[lambda] = avg_R2
    avg_train_correlation[lambda] = avg_correlation

    if (!is.null(Y_list_validation)) {

      R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, fit$Beta)
      correlation = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, fit$Beta)

      fit$performance$validation$R2 = R2
      fit$performance$validation$correlation = correlation

      error = compute_error(Y_list_validation, X_list_validation, indices_list_validation, as.matrix(fit$Beta))
      avg_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, fit$Beta)
      avg_correlation = compute_avg_correlation(Y_list_validation, X_list_validation, indices_list_validation, fit$Beta)

      validation_SSE[lambda] = error
      avg_validation_R2[lambda] = avg_R2
      avg_validation_correlation[lambda] = avg_correlation

    }

    models = c(models, list(fit))

  }

  tuning = list(train = list(SSE = train_SSE, avg_R2 = avg_train_R2, avg_correlation = avg_train_correlation),
                validation = list(SSE = validation_SSE, avg_R2 = avg_validation_R2, avg_correlation = avg_validation_correlation))

  return(list(model_fits = models,
              tuning = tuning,
              lambda_sequence = lambda_sequence,
              n_lambda = n_lambda))

}

#' @importFrom glmnet glmnet
fit_glmnet_univariate = function(Y, X, lambda_sequence) {

  fit = glmnet(X, Y, lambda = lambda_sequence/length(Y))

  return(as.matrix(coef(fit)))

}
