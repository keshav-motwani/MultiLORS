# evaluate_objective = function(Y_list, X_list, L_list, indices_list, Beta, lambda, gamma, gamma_weights) {
#
#   p1 = evaluate_g(Y_list, X_list, L_list, indices_list, Beta)
#
#   p2 = nuclear_norm_penalty(L_list, gamma, gamma_weights)
#
#   p3 = l1_penalty(Beta, lambda)
#
#   return(p1 + p2 + p3)
#
# }
#
# evaluate_g = function(Y_list, X_list, L_list, indices_list, Beta)  {
#
#   p1 = 0
#
#   for (k in 1:length(Y_list)) {
#
#     observed_indices = indices_list[[k]]
#
#     p1 = p1 + error(Y_list[[k]], (X_list[[k]] %*% Beta[, observed_indices]) + L_list[[k]])
#
#   }
#
#   return(p1 / 2)
#
# }
#
# l1_penalty = function(Beta, lambda) {
#
#   return(lambda * sum(abs(Beta[-1, ])))
#
# }
#
# nuclear_norm_penalty = function(L_list, gamma, gamma_weights) {
#
#   value = 0
#
#   for (k in 1:length(L_list)) {
#
#     gamma_k = gamma_weights[[k]] * gamma
#
#     decomp = svd(L_list[[k]])
#
#     value = value + sum(decomp$d) * gamma_k
#
#   }
#
#   return(value)
#
# }
#
# compute_error = function(Y_list, X_list, indices_list, Beta) {
#
#   error = 0
#
#   for (k in 1:length(Y_list)) {
#
#     error = error + sum((Y_list[[k]] - X_list[[k]] %*% Beta[, indices_list[[k]]]) ^ 2)
#
#   }
#
#   return(error)
#
# }

#' Compute SSE per response
#'
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Beta
#'
#' @return
#' @export
compute_SSE = function(Y_list, X_list, indices_list, Beta) {

    if (ncol(X_list[[1]]) + 1 == nrow(Beta)) {
      X_list = lapply(X_list, function(x) cbind(1, x))
    }

    error = numeric(ncol(Beta))

    for (k in 1:length(Y_list)) {

      error[indices_list[[k]]] = error[indices_list[[k]]] + colSums((Y_list[[k]] - X_list[[k]] %*% Beta[, indices_list[[k]]]) ^ 2)

    }

    return(error)

}

#' Compute SST per response
#'
#' @param Y_list_test
#' @param indices_list_test
#' @param Y_list_train
#' @param indices_list_train
#'
#' @return
#' @export
compute_SST = function(Y_list_test, indices_list_test, Y_list_train, indices_list_train) {

  q = max(unlist(indices_list_train))

  Y_train_means = numeric(q)
  counts = numeric(q)

  for (i in 1:q) {

    Y = subset_observed_data_univariate(Y_list_train, NULL, indices_list_train, i)$Y

    counts[i] = length(Y)

    Y_train_means[i] = mean(Y)

  }

  SST = numeric(q)

  for (k in 1:length(Y_list_test)) {

    preds = tcrossprod(rep(1, nrow(Y_list_test[[k]])), Y_train_means[indices_list_test[[k]]])

    SST[indices_list_test[[k]]] = SST[indices_list_test[[k]]] + colSums((Y_list_test[[k]] - preds) ^ 2)

  }

  attr(SST, "counts") = counts

  return(SST)

}

#' Compute R2 per response
#'
#' @param Y_list_test
#' @param X_list_test
#' @param indices_list_test
#' @param Y_list_train
#' @param indices_list_train
#' @param Beta
#'
#' @return
#' @export
compute_R2 = function(Y_list_test, X_list_test, indices_list_test, Y_list_train, indices_list_train, Beta) {

  SST = compute_SST(Y_list_test, indices_list_test, Y_list_train, indices_list_train)
  SSE = compute_SSE(Y_list_test, X_list_test, indices_list_test, Beta)

  return(1 - SSE/SST)

}

#' Compute average R2 across all responses
#'
#' @param Y_list_test
#' @param X_list_test
#' @param indices_list_test
#' @param Y_list_train
#' @param indices_list_train
#' @param Beta
#'
#' @return
#' @export
compute_avg_R2 = function(Y_list_test, X_list_test, indices_list_test, Y_list_train, indices_list_train, Beta) {

  SST = compute_SST(Y_list_test, indices_list_test, Y_list_train, indices_list_train)
  SSE = compute_SSE(Y_list_test, X_list_test, indices_list_test, Beta)

  return(mean(1 - SSE/SST, na.rm = TRUE))

}

#' Compute weighted average R2 across all responses with weights based on number of observations
#'
#' @param Y_list_test
#' @param X_list_test
#' @param indices_list_test
#' @param Y_list_train
#' @param indices_list_train
#' @param Beta
#'
#' @return
#' @export
compute_weighted_avg_R2 = function(Y_list_test, X_list_test, indices_list_test, Y_list_train, indices_list_train, Beta) {

  SST = compute_SST(Y_list_test, indices_list_test, Y_list_train, indices_list_train)
  SSE = compute_SSE(Y_list_test, X_list_test, indices_list_test, Beta)

  counts = attr(SST, "counts")
  counts[SST == 0] = 0
  weights = counts/sum(counts)

  return(sum(weights * (1 - SSE/SST), na.rm = TRUE))

}
