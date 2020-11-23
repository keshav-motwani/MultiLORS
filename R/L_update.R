update_L = function(Y_list, X_list, indices_list, Beta, gamma_weights, gamma) {

  L_list = vector(mode = "list", length = length(Y_list))

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    observed_indices = indices_list[[k]]

    L_list[[k]] = nuclear_prox(Y_list[[k]] - X_list[[k]] %*% Beta[, observed_indices], gamma_k)

  }

  return(L_list)

}