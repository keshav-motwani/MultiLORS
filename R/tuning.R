compute_gamma_weights = function(Y_list) {

  weights = sapply(Y_list, function(k) irlba::irlba(k, nv = 1)$d)

  return(weights)

}

compute_candidate_gamma_sequence = function(n_gamma, min_ratio) {

  gamma = log_seq(1, min_ratio, n_gamma)

  return(gamma)

}

compute_candidate_lambda_grid = function(Y_list, X_list, D_list, XtY_dot_list, gamma_weights, gamma_sequence, n_lambda, min_ratio) {

  lambda = sapply(gamma_sequence, function(gamma) compute_candidate_lambda_sequence_fixed_gamma(Y_list, X_list, D_list, XtY_dot_list, n_lambda, min_ratio, gamma_weights, gamma))

  return(lambda)

}

compute_candidate_lambda_sequence_fixed_gamma = function(Y_list, X_list, D_list, XtY_dot_list, n_lambda, min_ratio, gamma_weights, gamma) {

  result = matrix(0, nrow = nrow(XtY_dot_list[[1]]), ncol = ncol(XtY_dot_list[[1]]))

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    L_k = nuclear_prox(Y_list[[k]], gamma_k)

    result = result + crossprod(X_list[[k]], zero_pad_matrix(L_k, D_list[[k]])) - XtY_dot_list[[k]]

  }

  max_lambda = max(abs(result))

  lambda = log_seq(max_lambda, max_lambda * min_ratio, n_lambda)

  return(lambda)

}

log_seq = function(from, to, length) {

  exp(seq(log(from), log(to), length.out = length))

}
