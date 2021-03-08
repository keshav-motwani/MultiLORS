fit = function(Y_list,
               X_list,
               q,
               indices_list,
               XtX_list,
               XtY_list,
               lambda,
               gamma,
               gamma_weights,
               Beta_old,
               s_Beta,
               n_iter,
               tolerance,
               verbose) {

  objective = numeric(2 * n_iter)

  s = s_Beta * 100

  for (iter in 1:n_iter) {

    L_update = update_L(Y_list = Y_list,
                          X_list = X_list,
                          indices_list = indices_list,
                          Beta = Beta_old,
                          gamma_weights = gamma_weights,
                          gamma = gamma)
    L_list_new = L_update$L

    objective[2 * iter - 1] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_old, lambda, L_update$nuclear_norm_penalty)

    Beta_new = update_Beta(Y_list = Y_list,
                           X_list = X_list,
                           L_list = L_list_new,
                           q = q,
                           indices_list = indices_list,
                           XtX_list = XtX_list,
                           XtY_list = XtY_list,
                           Beta_old = Beta_old,
                           lambda = lambda,
                           s_Beta = s_Beta,
                           s = s)

    objective[2 * iter] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_new, lambda, L_update$nuclear_norm_penalty)

    last_step = objective[2 * (iter - 1)]
    this_step = objective[2 * iter]

    if (verbose == 2) print(paste0("Iteration ", iter, ": ", this_step))

    diff = (last_step - this_step)/last_step

    if (iter > 1 && diff < tolerance) {
      break
    }

    Beta_old = Beta_new

  }

  L_list_new = update_L(Y_list = Y_list,
                        X_list = X_list,
                        indices_list = indices_list,
                        Beta = Beta_new,
                        gamma_weights = gamma_weights,
                        gamma = gamma)$L

  if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda, " --- # of iterations: ", iter, "; difference = ", round((last_step - this_step)/last_step, 10)))

  objective = objective[1:(2 * iter)]

  result = list(
    Beta = Beta_new,
    L_list = L_list_new,
    objective = objective,
    n_iter = iter,
    lambda = lambda,
    gamma = gamma,
    diff = diff
  )

  return(result)

}

fit_solution_path = function(Y_list,
                             X_list,
                             indices_list,
                             Y_list_unstd,
                             Y_list_validation,
                             X_list_validation,
                             indices_list_validation,
                             standardize,
                             n_lambda,
                             n_gamma,
                             n_iter,
                             tolerance,
                             extra_iter,
                             extra_iter_threshold,
                             early_stopping,
                             verbose,
                             return_L,
                             p,
                             q,
                             XtX_list,
                             XtY_list,
                             X_mean,
                             X_sd,
                             Y_sd,
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
      verbose = verbose
    )

    model$lambda_index = lambda
    model$gamma_index = gamma

    Beta_old = model$Beta

    adjusted_Beta = adjust_Beta(model$Beta, X_mean, X_sd, Y_sd)
    colnames(adjusted_Beta) = attr(indices_list, "responses")
    if (!is.null(colnames(X_list[[1]]))) rownames(adjusted_Beta) = colnames(X_list[[1]])

    adjusted_L_list = adjust_L(model$L_list, indices_list, Y_sd)

    if (!is.null(Y_list_validation) & extra_iter > 0 & model$diff > tolerance) {

      avg_validation_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list_unstd, indices_list, adjusted_Beta)

      if (avg_validation_R2 > extra_iter_threshold) {

        extra_model = fit(
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
          n_iter = extra_iter,
          tolerance = tolerance,
          verbose = verbose
        )

        extra_model$objective = c(model$objective, extra_model$objective)
        extra_model$n_iter = model$n_iter + extra_model$n_iter

        extra_model$lambda_index = lambda
        extra_model$gamma_index = gamma

        model = extra_model

        Beta_old = model$Beta

        adjusted_Beta = adjust_Beta(model$Beta, X_mean, X_sd, Y_sd)
        colnames(adjusted_Beta) = attr(indices_list, "responses")
        if (!is.null(colnames(X_list[[1]]))) rownames(adjusted_Beta) = colnames(X_list[[1]])

        adjusted_L_list = adjust_L(model$L_list, indices_list, Y_sd)

      }

    }

    model$performance = list(train = list(), validation = list())

    model$performance$train$R2 = compute_R2(Y_list, X_list, indices_list, Y_list, indices_list, model$Beta)
    model$performance$train$correlation = compute_correlation(Y_list, X_list, indices_list, model$Beta)

    if (!is.null(Y_list_validation)) {

      validation_error = compute_error(Y_list_validation, X_list_validation, indices_list_validation, adjusted_Beta)
      avg_validation_R2 = compute_avg_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list_unstd, indices_list, adjusted_Beta)

      min_validation_error = min(min_validation_error, validation_error)
      max_avg_validation_R2 = max(max_avg_validation_R2, avg_validation_R2)

      model$performance$validation$R2 = compute_R2(Y_list_validation, X_list_validation, indices_list_validation, Y_list_unstd, indices_list, adjusted_Beta)
      model$performance$validation$correlation  = compute_correlation(Y_list_validation, X_list_validation, indices_list_validation, adjusted_Beta)

      if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda, " --- Validation Error: ", validation_error, "; Avg Validation R2: ", avg_validation_R2))

    }

    model$Beta = as(adjusted_Beta, "dgCMatrix")

    if (return_L) {
      model$L_list = adjusted_L_list
    } else {
      model$L_list = NULL
    }

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