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
                     indices_list,
                     Y_list_validation = NULL,
                     X_list_validation = NULL,
                     indices_list_validation = NULL,
                     n_lambda = 20,
                     n_gamma = 20,
                     lambda_min_ratio = 0.001,
                     gamma_min_ratio = 0.001,
                     n_iter = 1000,
                     tolerance = 1e-6,
                     early_stopping = TRUE,
                     verbose = 0,
                     return_L = TRUE,
                     n_cores = 1) {

  X_list_unstd = X_list
  X_list = standardize_X(X_list)
  X_mean = attributes(X_list)$mean
  X_sd = attributes(X_list)$sd
  X_list = lapply(X_list, function(k) cbind(1, k))
  XtX_list = lapply(X_list, crossprod)
  XtY_list = mapply(crossprod, x = X_list, y = Y_list, SIMPLIFY = FALSE)
  p = ncol(X_list[[1]])
  q = max(unlist(indices_list))
  K = length(X_list)
  dataset_indices_list = lapply(1:q, function(j) which(sapply(1:K, function(k) j %in% indices_list[[k]])))
  if (!is.null(X_list_validation)) {
    X_list_validation = lapply(X_list_validation, function(k) cbind(1, k))
  }

  gamma_weights = compute_gamma_weights(Y_list)
  gamma_sequence = compute_candidate_gamma_sequence(n_gamma, gamma_min_ratio)
  lambda_grid = compute_candidate_lambda_grid(Y_list, X_list, q, indices_list, XtY_list, gamma_weights, gamma_sequence, n_lambda, lambda_min_ratio)

  s_Beta = compute_s_Beta(XtX_list, p, q, dataset_indices_list)

  model_fits = parallel::mclapply(1:n_gamma, function(gamma) {
    fit_solution_path(
      Y_list,
      X_list,
      indices_list,
      Y_list_validation,
      X_list_validation,
      indices_list_validation,
      standardize,
      n_lambda,
      n_gamma,
      lambda_min_ratio,
      gamma_min_ratio,
      n_iter,
      tolerance,
      early_stopping,
      verbose,
      return_L,
      p,
      q,
      XtX_list,
      XtY_list,
      X_mean,
      X_sd,
      lambda_grid,
      gamma_sequence,
      gamma_weights,
      s_Beta,
      gamma
    )
  }, mc.cores = n_cores
  )

  fit = list(model_fits = model_fits,
             n_lambda = n_lambda,
             n_gamma = n_gamma,
             gamma_sequence = gamma_sequence,
             gamma_weights = gamma_weights,
             lambda_grid = lambda_grid)

  train = compute_tuning_performance(fit, Y_list, X_list_unstd, indices_list, Y_list, indices_list)

  if (!is.null(X_list_validation)) {
    validation = compute_tuning_performance(fit, Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list)
  } else {
    validation = NULL
  }

  fit = c(fit, list(tuning = list(train = train, validation = validation)))

  return(fit)

}

#' Refit MultiLORS on nonzero coefficients and corrected responses using OLS
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
refit_MultiLORS = function(fit,
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

  refit_Betas = parallel::mclapply(
    fit$model_fits,
    function(solution_path) {
      print(paste0("gamma: ", solution_path[[1]]$gamma_index))
      lapply(solution_path, function(model) {
        print(paste0("gamma: ", model$gamma_index, "lambda: ", model$lambda_index))
        MultiLORS:::refit_OLS(Y_list, X_list, model$L_list, indices_list, model$Beta)
      })
    },
  mc.cores = n_cores)

  for (solution_path in fit$model_fits) {

    for (model in solution_path) {

      gamma = model$gamma_index
      lambda = model$lambda_index
      L_list = model$L_list
      Beta = model$Beta

      refit_Beta = refit_Betas[[gamma]][[lambda]]
      fit$model_fits[[gamma]][[lambda]]$Beta = refit_Beta
      refit_Beta = as.matrix(refit_Beta)

      fit$model_fits[[gamma]][[lambda]]$performance$train$R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, refit_Beta)
      fit$model_fits[[gamma]][[lambda]]$performance$train$correlation = compute_correlation(Y_list, X_list, indices_list, refit_Beta)

      if (!is.null(Y_list_validation)) {
        fit$model_fits[[gamma]][[lambda]]$performance$validation$R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, refit_Beta)
        fit$model_fits[[gamma]][[lambda]]$performance$validation$correlation  = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, refit_Beta)
      }

    }

  }

  train = compute_tuning_performance(fit, Y_list, X_list, indices_list, Y_list, indices_list)

  if (!is.null(X_list_validation)) {
    validation = compute_tuning_performance(fit, Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list)
  } else {
    validation = NULL
  }

  fit$tuning = list(train = train, validation = validation)

  return(fit)

}
