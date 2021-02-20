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

  models = list()

  for (lambda in 1:n_lambda) {

    fit = list(Beta = as(Beta[, , lambda], "dgCMatrix"), performance = list(train = list(), validation = list()))
    colnames(fit$Beta) = attr(indices_list, "responses")
    if (!is.null(colnames(X_list[[1]]))) rownames(fit$Beta) = colnames(X_list[[1]])

    R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, fit$Beta)
    correlation = compute_correlation(Y_list, X_list, indices_list, fit$Beta)

    fit$performance$train$R2 = R2
    fit$performance$train$correlation = correlation

    if (!is.null(Y_list_validation)) {

      R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, fit$Beta)
      correlation = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, fit$Beta)

      fit$performance$validation$R2 = R2
      fit$performance$validation$correlation = correlation

    }

    models = c(models, list(fit))

  }

  fit = list(model_fits = models,
             lambda_sequence = lambda_sequence,
             n_lambda = n_lambda)

  train = compute_tuning_performance_glmnet(fit, Y_list, X_list, indices_list, Y_list, indices_list)

  if (!is.null(X_list_validation)) {
    validation = compute_tuning_performance_glmnet(fit, Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list)
  } else {
    validation = NULL
  }

  fit$tuning = list(train = train, validation = validation)

  return(fit)

}

#' Refit glmnet on nonzero coefficients using OLS
#'
#' @param fit
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Y_list_validation
#' @param X_list_validation
#' @param indices_list_validation
#'
#' @return
#' @export
refit_glmnet = function(fit,
                        Y_list,
                        X_list,
                        indices_list,
                        Y_list_validation = NULL,
                        X_list_validation = NULL,
                        indices_list_validation = NULL,
                        n_cores = 1) {

  X_list = lapply(X_list, function(k) cbind(1, k))

  if (!is.null(X_list_validation)) {
    X_list_validation = lapply(X_list_validation, function(k) cbind(1, k))
  }

  q = max(unlist(indices_list))

  subsetted_X = lapply(1:q, function(i) subset_observed_data_univariate(NULL, X_list, indices_list, i)$X)
  subsetted_XtX = lapply(subsetted_X, crossprod)

  refit_Betas = parallel::mclapply(
    fit$model_fits,
    function(model) {
      refit_OLS(Y_list, subsetted_XtX, subsetted_X, NULL, indices_list, model$Beta)
    },
    mc.cores = n_cores)

  for (lambda in 1:length(fit$model_fits)) {

    Beta = fit$model_fits[[lambda]]$Beta

    refit_Beta = refit_Betas[[lambda]]
    fit$model_fits[[lambda]]$Beta = refit_Beta
    refit_Beta = as.matrix(refit_Beta)

    fit$model_fits[[lambda]]$performance$train$R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, refit_Beta)
    fit$model_fits[[lambda]]$performance$train$correlation = compute_correlation(Y_list, X_list, indices_list, refit_Beta)

    if (!is.null(Y_list_validation)) {
      fit$model_fits[[lambda]]$performance$validation$R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, refit_Beta)
      fit$model_fits[[lambda]]$performance$validation$correlation  = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, refit_Beta)
    }

  }

  train = compute_tuning_performance_glmnet(fit, Y_list, X_list, indices_list, Y_list, indices_list)

  if (!is.null(X_list_validation)) {
    validation = compute_tuning_performance_glmnet(fit, Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list)
  } else {
    validation = NULL
  }

  fit$tuning = list(train = train, validation = validation)

  return(fit)

}

#' Compute tuning performance on training/validation data from glmnet fit
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
compute_tuning_performance_glmnet = function(fit, Y_list, X_list, indices_list, Y_list_train, indices_list_train) {

  n_lambda = fit$n_lambda

  p = nrow(fit[[1]]$Beta)
  q = ncol(fit[[1]]$Beta)

  if (ncol(X_list[[1]]) + 1 == nrow(fit$model_fits[[1]]$Beta)) {
    X_list = lapply(X_list, function(x) cbind(1, x))
  }

  SSE = numeric(n_lambda)
  avg_R2 = numeric(n_lambda)
  avg_correlation = numeric(n_lambda)

  for (i in 1:length(fit$model_fits)) {

      Beta = as.matrix(fit$model_fits[[i]]$Beta)

      SSE[i] = compute_error(Y_list, X_list, indices_list, Beta)
      avg_R2[i] = compute_avg_R2(Y_list, X_list, indices_list, Y_list_train, indices_list_train, Beta)
      avg_correlation[i] = compute_avg_correlation(Y_list, X_list, indices_list, Beta)

  }

  return(list(SSE = SSE,
              avg_R2 = avg_R2,
              avg_correlation = avg_correlation))

}

#' @importFrom glmnet glmnet
fit_glmnet_univariate = function(Y, X, lambda_sequence) {

  fit = glmnet(X, Y, lambda = lambda_sequence/length(Y))

  return(as.matrix(coef(fit)))

}
