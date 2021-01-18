update_L = function(Y_list, X_list, indices_list, Beta, gamma_weights, gamma) {

  L_list = vector(mode = "list", length = length(Y_list))

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    observed_indices = indices_list[[k]]

    L_list[[k]] = nuclear_prox(Y_list[[k]] - X_list[[k]] %*% Beta[, observed_indices], gamma_k)

  }

  return(L_list)

}

compress_L = function(L_list) {

  result = list()

  for (i in 1:length(L_list)) {

    decomp = svd(L_list[[i]])

    indices = decomp$d > 1e-10

    decomp = list(u = decomp$u[, indices, drop = FALSE],
                  d = decomp$d[indices],
                  v = decomp$v[, indices, drop = FALSE])

    result = c(result, list(decomp))

  }

  return(result)

}