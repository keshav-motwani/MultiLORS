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
               s_Beta,
               n_iter,
               tolerance,
               verbose,
               return_L) {

  objective = numeric(2 * n_iter)

  if (is.null(Beta_old)) Beta_old = matrix(0, nrow = ncol(X_list[[1]]), ncol = q)

  s = s_Beta * 100

  for (iter in 1:n_iter) {

    L_list_new = update_L(Y_list = Y_list,
                          X_list = X_list,
                          indices_list = indices_list,
                          Beta = Beta_old,
                          gamma_weights = gamma_weights,
                          gamma = gamma)
    objective[2 * iter - 1] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_old, lambda, gamma, gamma_weights)

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
    objective[2 * iter] = evaluate_objective(Y_list, X_list, L_list_new, indices_list, Beta_new, lambda, gamma, gamma_weights)

    last_step = objective[2 * (iter - 1)]
    this_step = objective[2 * iter]

    if (verbose == 2) print(paste0("Iteration ", iter, ": ", this_step))

    if (iter > 1 && (last_step - this_step)/last_step < tolerance) {
      break
    }

    Beta_old = Beta_new
    L_list_old = L_list_new

  }

  if (verbose > 0) print(paste0("gamma: ", gamma, "; lambda: ", lambda, " --- # of iterations: ", iter, "; difference = ", round((last_step - this_step)/last_step, 5)))

  objective = objective[1:(2 * iter)]

  if (return_L) {
    L_list = compress_L(L_list_new)
  } else {
    L_list = NULL
  }

  result = list(
    Beta = adjust_Beta(Beta_new, X_mean, X_sd),
    L_list = L_list,
    objective = objective,
    n_iter = iter,
    lambda = lambda,
    gamma = gamma
  )

  return(result)

}