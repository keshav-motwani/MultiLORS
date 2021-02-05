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

  if (!is.matrix(Beta)) {
    Beta = as.matrix(Beta)
  }

  error = numeric(ncol(Beta))

  for (k in 1:length(Y_list)) {

    error[indices_list[[k]]] = error[indices_list[[k]]] + colSums((Y_list[[k]] - X_list[[k]] %*% Beta[, indices_list[[k]]]) ^ 2)

  }

  names(error) = colnames(Beta)

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

#' Compute correlation between actual and predicted per response
#'
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Beta
#'
#' @return
#' @export
compute_correlation = function(Y_list, X_list, indices_list, Beta) {

  if (ncol(X_list[[1]]) + 1 == nrow(Beta)) {
    X_list = lapply(X_list, function(x) cbind(1, x))
  }

  if (!is.matrix(Beta)) {
    Beta = as.matrix(Beta)
  }

  corr = matrix(nrow = length(Y_list), ncol = ncol(Beta))

  for (k in 1:length(Y_list)) {

    pred = X_list[[k]] %*% Beta[, indices_list[[k]]]

    for (i in indices_list[[k]]) {

      correlation = cor(Y_list[[k]][, i], pred[, i])

      corr[k, i] = ifelse(is.na(correlation), 0, correlation)

    }

  }

  result = colMeans(corr, na.rm = TRUE)

  names(result) = colnames(Beta)

  return(result)

}

#' Compute average correlation across all responses
#'
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Beta
#'
#' @return
#' @export
compute_avg_correlation = function(Y_list, X_list, indices_list, Beta) {

  corr = compute_correlation(Y_list, X_list, indices_list, Beta)

  return(mean(corr, na.rm = TRUE))

}

#' Determine indices of minimum element in matrix
#'
#' @param mat
#'
#' @return
#' @export
which_min = function(mat) {

  which(mat == min(mat, na.rm = TRUE), arr.ind = TRUE)

}

#' Determine indices of maximum element in matrix
#'
#' @param mat
#'
#' @return
#' @export
which_max = function(mat) {

  which(mat == max(mat, na.rm = TRUE), arr.ind = TRUE)

}