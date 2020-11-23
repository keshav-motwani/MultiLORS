update_Beta = function(Y_list, X_list, L_list, D_list, indices_list, XtX_list, XtY_dot_list, Beta_old, lambda, s_Beta, s, line_search) {

  gradient = compute_gradient_Beta(X_list = X_list,
                                   L_list = L_list,
                                   D_list = D_list,
                                   XtX_list = XtX_list,
                                   XtY_dot_list = XtY_dot_list,
                                   Beta_old = Beta_old)

  shrinkage = 0.5

  # https://people.eecs.berkeley.edu/~elghaoui/Teaching/EE227A/lecture18.pdf

  if (line_search) {

    g_old = evaluate_g(Y_list, X_list, L_list, indices_list, Beta_old)

    while(line_search & s != s_Beta) {

      Beta = l1_prox(Beta_old - (s * gradient), s * lambda)

      g_new = evaluate_g(Y_list, X_list, L_list, indices_list, Beta)

      difference = (Beta_old - Beta) / s

      if (g_new > g_old - s * sum(gradient * difference) + 0.5 * s * sum(difference^2)) {

        s = max(shrinkage * s, s_Beta)

      } else {

        line_search = FALSE

      }

    }

  } else {

    Beta = l1_prox(Beta_old - (s_Beta * gradient), s_Beta * lambda)
    s = 0

  }

  # print(paste0("s: ", s, "; s_Beta: ", s_Beta))

  return(list(Beta = Beta, s = s))

}

compute_gradient_Beta = function(X_list, L_list, D_list, XtX_list, XtY_dot_list, Beta_old) {

  gradient = matrix(0, nrow = ncol(X_list[[1]]), ncol = nrow(D_list[[1]]))

  for (k in 1:length(XtY_dot_list)) {

    one = XtX_list[[k]] %*% Beta_old %*% D_list[[k]]

    two = crossprod(X_list[[k]], zero_pad_matrix(L_list[[k]], D_list[[k]]))

    three = XtY_dot_list[[k]]

    gradient = gradient + one + two - three

  }

  return(gradient)

}

compute_s_Beta = function(XtX_list, D_list) {

  max_eigenvalues = numeric(nrow(D_list[[1]]))

  for (j in 1:length(max_eigenvalues)) {

    P_j = which(sapply(D_list, function(k) k[j, j] == 1))

    block = Reduce('+', XtX_list[P_j])

    max_eigenvalues[j] = eigen(block, TRUE)$values[1]

  }

  return(1 / max(max_eigenvalues))

}
