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
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
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

  # https://github.com/r-pkg-examples/rcpp-and-doparallel
  cl = makeCluster(n_cores, outfile="")
  on.exit(stopCluster(cl))
  registerDoParallel(cl)

  result = foreach(gamma = 1:n_gamma) %dopar% {
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
  }

  result = aggregate_tuning_results(result, n_lambda, n_gamma)

  result$gamma_sequence = gamma_sequence
  result$lambda_grid = lambda_grid

  return(result)

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

  for (lambda in 1:n_lambda) {

    if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda))

    model = fit(
      Y_list = Y_list,
      X_list = X_list,
      q = q,
      indices_list = indices_list,
      XtX_list = XtX_list,
      XtY_list = XtY_list,
      X_mean = X_mean,
      X_sd = X_sd,
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

    if (!is.null(Y_list_validation)) {
      model$validation_error = compute_error(Y_list_validation, X_list_validation, indices_list_validation, model$Beta)
      model$avg_validation_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, model$Beta)
      model$weighted_avg_validation_R2 = compute_weighted_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, model$Beta)
      model$R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list, indices_list, model$Beta)
      if (verbose > 0) print(paste0("Validation Error: ", model$validation_error, "; Avg Validation R2: ", round(model$avg_validation_R2, 4), "; Weighted Avg Validation R2: ", round(model$weighted_avg_validation_R2, 4)))
    }

    Beta_old = model$Beta
    L_list_old = model$L_list

    model$Beta = as(model$Beta, "dgCMatrix")
    colnames(model$Beta) = attr(indices_list, "responses")

    result = c(result, list(model))

    if (early_stopping && !is.null(Y_list_validation)) {
      if (lambda > 5 &&
          lambda > n_lambda / 4 &&
          model$validation_error > result[[length(result) - 1]][["validation_error"]] &&
          model$avg_validation_R2 < result[[length(result) - 1]][["avg_validation_R2"]] &&
          model$weighted_avg_validation_R2 < result[[length(result) - 1]][["weighted_avg_validation_R2"]]) {
        break
      }
    }


  }

  return(result)

}
