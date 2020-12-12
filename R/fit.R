fit = function(Y_list,
               X_list,
               q,
               indices_list,
               XtX_list,
               XtY_list,
               X_mean,
               X_sd,
               lambda,
               gamma,
               gamma_weights,
               Beta_old,
               L_list_old,
               s_Beta,
               n_iter,
               tolerance,
               line_search,
               verbose) {

  objective = numeric(2 * n_iter)

  if (is.null(Beta_old)) Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = q)
  if (is.null(L_list_old)) L_list_old = lapply(Y_list, function(k) matrix(0, nrow = nrow(k), ncol = ncol(k)))

  s = s_Beta * 10

  for (iter in 1:n_iter) {

    L_list_new = update_L(Y_list = Y_list,
                          X_list = X_list,
                          indices_list = indices_list,
                          Beta = Beta_old,
                          gamma_weights = gamma_weights,
                          gamma = gamma)
    objective[2 * iter - 1] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_old, lambda, gamma, gamma_weights)

    Beta_update = update_Beta(Y_list = Y_list,
                              X_list = X_list,
                              L_list = L_list_new,
                              q = q,
                              indices_list = indices_list,
                              XtX_list = XtX_list,
                              XtY_list = XtY_list,
                              Beta_old = Beta_old,
                              lambda = lambda,
                              s_Beta = s_Beta,
                              s = s,
                              line_search = line_search)
    Beta_new = Beta_update$Beta
    objective[2 * iter] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_new, lambda, gamma, gamma_weights)

    last_step = objective[2 * (iter - 1)]
    this_step = objective[2 * iter]

    if (verbose) print(paste0("Iteration ", iter, ": ", this_step))

    if (iter > 1 && (last_step - this_step)/last_step < tolerance) {
      break
    }

    Beta_old = Beta_new
    L_list_old = L_list_new
    s = Beta_update$s * 10

  }

  objective = objective[1:(2 * iter)]

  result = list(
    Beta = adjust_Beta(Beta_new, X_mean, X_sd),
    L_list = L_list_new,
    objective = objective,
    n_iter = iter,
    lambda = lambda,
    gamma = gamma
  )

  return(result)

}