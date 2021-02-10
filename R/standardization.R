standardize_X = function(X_list) {

  X = do.call(rbind, X_list)
  k = rep(1:length(X_list), sapply(X_list, nrow))

  means = Matrix::colMeans(X)
  vars = matrixStats::colVars(as.matrix(X)) * (nrow(X) - 1) / nrow(X)

  X = X - tcrossprod(rep(1, nrow(X)), means)
  X = X %*% diag(1 / sqrt(vars))
  colnames(X) = colnames(X_list[[1]])

  X_list = lapply(1:length(X_list), function(i) X[k == i, ])

  attributes(X_list)$mean = means
  attributes(X_list)$sd = sqrt(vars)

  return(X_list)

}

adjust_Beta = function(Beta, X_mean, X_sd) {

  if (!is.null(X_mean)) {

    Beta[1, ] = Beta[1, ] - crossprod(X_mean / X_sd, Beta[-1, ])

    Beta[-1, ] = diag(1 / X_sd) %*% Beta[-1, ]

  }

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
  # Beta1[1, ] = Beta1[1, ] - v
  #
  # Beta1[-1, ] = diag(1 / X_sd) %*% Beta1[-1, ]
  # head(Beta1)
  #
  # model = glmnet(X, Y, family = "mgaussian", standardize = TRUE)
  # Beta2 = aperm(simplify2array(lapply(coef(model), as.matrix)), c(2, 1, 3))[3, , ]
  # head(Beta2)

}
