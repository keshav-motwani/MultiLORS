#' @importFrom matrixStats colVars
standardize_X = function(X_list) {

  X = do.call(rbind, X_list)
  k = rep(1:length(X_list), sapply(X_list, nrow))

  means = colMeans(X)
  vars = colVars(as.matrix(X)) * (nrow(X) - 1) / nrow(X)

  X = X - tcrossprod(rep(1, nrow(X)), means)
  X = X %*% diag(1 / sqrt(vars))
  colnames(X) = colnames(X_list[[1]])

  X_list = lapply(1:length(X_list), function(i) X[k == i, ])

  attributes(X_list)$mean = means
  attributes(X_list)$sd = sqrt(vars)

  return(X_list)

}

#' @importFrom matrixStats colVars
standardize_Y = function(Y_list, indices_list, dataset_indices_list) {

  q = max(unlist(indices_list))

  sds = lapply(Y_list, function(Y) colVars(Y))

  n = lapply(Y_list, nrow)

  pooled_sds = numeric(q)
  pooled_n = numeric(q)
  n_datasets = numeric(q)

  for (k in 1:length(Y_list)) {

    pooled_sds[indices_list[[k]]] = pooled_sds[indices_list[[k]]] + sds[[k]] * (n[[k]] - 1)
    pooled_n[indices_list[[k]]] = pooled_n[indices_list[[k]]] + n[[k]]
    n_datasets[indices_list[[k]]] = n_datasets[indices_list[[k]]] + 1

  }

  pooled_sds = sqrt(pooled_sds / (pooled_n - n_datasets))

  Y_list = lapply(Y_list, function(Y) Y %*% diag(x = 1/pooled_sds[indices_list[[k]]], nrow = ncol(Y)))

  attributes(Y_list)$sd = pooled_sds

  return(Y_list)

}

adjust_Beta = function(Beta, X_mean, X_sd, Y_sd) {

  Beta[1, ] = Beta[1, ] - crossprod(X_mean / X_sd, Beta[-1, ])

  Beta[-1, ] = diag(1 / X_sd) %*% Beta[-1, ]

  Beta = Beta %*% diag(x = Y_sd, nrow = length(Y_sd))

  return(Beta)

  # X = matrix((rnorm(100 * 20) + 3) / 7, 100, 20)
  # Y = matrix(rnorm(100 * 3), 100, 3)
  #
  # X_centered <- apply(X, 2, function(x) x - mean(x))
  # Xs <- apply(X_centered, 2, function(x) x / sqrt(sum(x^2) / nrow(X)))
  #
  # X_mean <- colMeans(X)
  # X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / nrow(X)))
  #
  # model = glmnet(Xs, Y, family = "mgaussian", standardize = FALSE)
  # Beta1 = aperm(simplify2array(lapply(coef(model), as.matrix)), c(2, 1, 3))[3, , ]
  #
  # Beta1[1, ] = Beta1[1, ] - crossprod(X_mean / X_sd, Beta1[-1, ])
  #
  # Beta1[-1, ] = diag(1 / X_sd) %*% Beta1[-1, ]
  # head(Beta1)
  #
  # model = glmnet(X, Y, family = "mgaussian", standardize = TRUE)
  # Beta2 = aperm(simplify2array(lapply(coef(model), as.matrix)), c(2, 1, 3))[3, , ]
  # head(Beta2)

}

adjust_L = function(L_list, indices_list, Y_sd) {

  mapply(L_list, indices_list, FUN = function(L, indices) L %*% diag(x = Y_sd[indices], nrow = ncol(L)), SIMPLIFY = FALSE)

}
