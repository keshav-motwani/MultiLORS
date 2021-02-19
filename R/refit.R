refit_OLS = function(Y_list, subsetted_X, L_list, indices_list, Beta) {

  if (!is.null(L_list)) {
    Y_list = mapply(Y = Y_list, L = L_list, FUN = function(Y, L) {
      D = diag(length(L$d))
      diag(D) = L$d
      Y - L$u %*% D %*% t(L$v)
    }, SIMPLIFY = FALSE)
  }

  p = nrow(Beta)
  q = ncol(Beta)

  refitted_Beta = matrix(nrow = p, ncol = q)

  for (i in 1:q) {

    # print(i)

    coefficients = numeric(p)

    nonzero = unique(c(1, which(Beta[, i] != 0)))

    data = subset_observed_data_univariate(Y_list, NULL, indices_list, i)

    coefficients[nonzero] = OLS(subsetted_X[[i]][, nonzero, drop = FALSE], data$Y)

    refitted_Beta[, i] = coefficients

  }

  colnames(refitted_Beta) = colnames(Beta)
  rownames(refitted_Beta) = rownames(Beta)

  return(as(refitted_Beta, "dgCMatrix"))

}
