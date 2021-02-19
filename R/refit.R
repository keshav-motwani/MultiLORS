refit_OLS = function(Y_list, X_list, L_list, indices_list, Beta) {

  if (!is.null(L_list)) {
    Y_list = mapply(Y = Y_list, L = L_list, FUN = function(Y, L) {
      D = diag(length(L$d))
      diag(D) = L$d
      Y - L$u %*% D %*% t(L$v)
    }, SIMPLIFY = FALSE)
  }

  p = ncol(X_list[[1]])
  q = max(unlist(indices_list))

  refitted_Beta = matrix(nrow = p, ncol = q)

  for (i in 1:q) {

    print(i)

    coefficients = numeric(p)

    nonzero = unique(c(1, which(Beta[, i] != 0)))

    data = subset_observed_data_univariate(Y_list, X_list, indices_list, i)

    coefficients[nonzero] = OLS(data$X[, nonzero, drop = FALSE], data$Y)

    refitted_Beta[, i] = coefficients

  }

  colnames(refitted_Beta) = colnames(Beta)
  rownames(refitted_Beta) = rownames(Beta)

  return(as(refitted_Beta, "dgCMatrix"))

}
