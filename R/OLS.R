#' Fit OLS on concatenated data
#'
#' @param Y_list
#' @param X_list
#' @param indices_list
#'
#' @return
#' @export
fit_OLS = function(Y_list,
                   X_list,
                   indices_list) {

  p = ncol(X_list[[1]])
  q = max(unlist(indices_list))

  Beta = matrix(nrow = p + 1, ncol = q)

  for (index in 1:q) {

    data = subset_observed_data_univariate(Y_list, X_list, indices_list, index)

    data$X = cbind(1, data$X)

    Beta[, index] = OLS(crossprod(data$X), data$X, data$Y)

  }

  return(Beta)

}