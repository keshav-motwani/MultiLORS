#' Construct D_{(k)} matrices based on observed responses
#'
#' @param observed_response_list
#'
#' @return
#' @export
construct_D_list = function(Y_list) {

  union = sort(unique(unlist(lapply(Y_list, colnames))))

  map = 1:length(union)
  names(map) = union

  for (k in 1:length(Y_list)) {

    Y_list[[k]] = Y_list[[k]][, order(colnames(Y_list[[k]])), drop = FALSE]

  }

  indices = lapply(Y_list, function(k) map[colnames(k)])

  D_list = lapply(indices, function(i) construct_D_from_indices(i, length(union)))

  for (k in 1:length(D_list)) {

    colnames(D_list[[k]]) = union
    rownames(D_list[[k]]) = union

  }

  return(list(Y_list = Y_list, D_list = D_list))

}

construct_D_from_indices = function(observed_indices, q) {

  D = matrix(0, nrow = q, ncol = q)
  diag(D)[observed_indices] = 1

  return(D)

}

zero_pad_matrix = function(matrix, D) {

  q = nrow(D)

  observed_indices = which(diag(D) == 1)

  padded = matrix(0, nrow = nrow(matrix), ncol = q)
  padded[, observed_indices] = matrix

  return(padded)

}

subset_observed_data_univariate = function(Y_list, X_list, D_list, index) {

  new_Y_list = list()
  new_X_list = list()

  for (k in 1:length(Y_list)) {

    if (D_list[[k]][index, index] == 1) {

      new_Y_list = c(new_Y_list, list(zero_pad_matrix(Y_list[[k]], D_list[[k]])[, index]))
      new_X_list = c(new_X_list, list(X_list[[k]]))

    }

  }

  Y = do.call(c, new_Y_list)
  X = do.call(rbind, new_X_list)
  batch = rep(1:length(new_Y_list), sapply(new_Y_list, length))

  return(list(Y = Y, X = X, batch = batch))

}
