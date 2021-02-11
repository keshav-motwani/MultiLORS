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

fit_solution_path = function(Y_list,
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
                             gamma) {

  result = list()

  Beta_old = matrix(0, nrow = p, ncol = q)

  min_validation_error = Inf
  max_avg_validation_R2 = 0

  for (lambda in 1:n_lambda) {

    if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda))

    model = fit(
      Y_list = Y_list,
      X_list = X_list,
      q = q,
      indices_list = indices_list,
      XtX_list = XtX_list,
      XtY_list = XtY_list,
      lambda = lambda_grid[gamma, lambda],
      gamma = gamma_sequence[gamma],
      gamma_weights = gamma_weights,
      Beta_old = Beta_old,
      s_Beta = s_Beta,
      n_iter = n_iter,
      tolerance = tolerance,
      verbose = verbose,
      return_L
    )

    model$lambda_index = lambda
    model$gamma_index = gamma

    model$performance = list(train = list(), validation = list())

    model$performance$train$R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, model$Beta)
    model$performance$train$correlation = compute_correlation(Y_list, X_list, indices_list, model$Beta)

    adjusted_Beta = adjust_Beta(model$Beta, X_mean, X_sd)

    if (!is.null(Y_list_validation)) {
      validation_error = compute_error(Y_list_validation, X_list_validation, indices_list_validation, adjusted_Beta)
      avg_validation_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, adjusted_Beta)
      min_validation_error = min(min_validation_error, validation_error)
      max_avg_validation_R2 = max(max_avg_validation_R2, avg_validation_R2)
      model$performance$validation$R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, adjusted_Beta)
      model$performance$validation$correlation  = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, adjusted_Beta)
      if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda, " --- Validation Error: ", validation_error, "; Avg Validation R2: ", round(avg_validation_R2, 4)))
    }

    Beta_old = model$Beta
    L_list_old = model$L_list

    model$Beta = as(adjusted_Beta, "dgCMatrix")
    colnames(model$Beta) = attr(indices_list, "responses")
    if (!is.null(colnames(X_list[[1]]))) rownames(model$Beta) = colnames(X_list[[1]])

    result = c(result, list(model))

    if (early_stopping && !is.null(Y_list_validation)) {
      if (lambda > 5 &&
          lambda > n_lambda / 4 &&
          validation_error > min_validation_error * 1.01 &&
          avg_validation_R2 < max_avg_validation_R2 * 0.99) {
        break
      }
    }

  }

  return(result)

}
