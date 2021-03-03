update_L = function(Y_list, X_list, indices_list, Beta, gamma_weights, gamma) {

  L_list = vector(mode = "list", length = length(Y_list))
  nuclear_norm_penalty = 0

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    observed_indices = indices_list[[k]]

    prox = nuclear_prox(Y_list[[k]] - X_list[[k]] %*% Beta[, observed_indices], gamma_k)

    L_list[[k]] = prox$L

    nuclear_norm_penalty = nuclear_norm_penalty + prox$nuclear_norm_penalty

  }

  return(list(L = L_list, nuclear_norm_penalty = nuclear_norm_penalty))

}

#' Recompute L
#'
#' @param fit
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param gamma
#' @param lambda
#'
#' @return
#' @export
recompute_L = function(fit, Y_list, X_list, indices_list, gamma, lambda) {

  Beta = as.matrix(fit$model_fits[[gamma]][[lambda]]$Beta)
  gamma_weights = fit$gamma_weights
  gamma = fit$gamma_sequence[gamma]
  Y_sd = fit$standardization$Y_sd

  if (ncol(X_list[[1]]) + 1 == nrow(Beta)) {
    X_list = lapply(X_list, function(x) cbind(1, x))
  }

  L_list = vector(mode = "list", length = length(Y_list))

  for (k in 1:length(Y_list)) {

    gamma_k = gamma_weights[[k]] * gamma

    observed_indices = indices_list[[k]]

    prox = nuclear_prox((Y_list[[k]] - X_list[[k]] %*% Beta[, observed_indices]) %*% diag(1 / Y_sd), gamma_k)

    L_list[[k]] = prox$L %*% diag(Y_sd)

  }

  return(L_list)

}