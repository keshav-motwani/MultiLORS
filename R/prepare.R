#' Make ordering of observed responses consistent and construct vectors indicating indices of all responses observed
#'
#' @param Y_list_train
#'
#' @return
#' @export
prepare_Y_and_indices_train = function(Y_list_train) {

  union = sort(unique(unlist(lapply(Y_list_train, colnames))))

  map = 1:length(union)
  names(map) = union

  for (k in 1:length(Y_list_train)) {

    Y_list_train[[k]] = Y_list_train[[k]][, sort(colnames(Y_list_train[[k]])), drop = FALSE]

  }

  indices = lapply(Y_list_train, function(k) map[colnames(k)])

  attr(indices, "map") = map
  attr(indices, "responses") = names(map)

  return(list(Y_list = Y_list_train, indices_list = indices))

}

#' Make ordering of observed responses consistent and construct vectors indicating indices of all responses observed, based on responses available in training data
#'
#' @param Y_list_test
#'
#' @return
#' @export
prepare_Y_and_indices_test = function(Y_list_test, indices_list_train) {

  map = attr(indices_list_train, "map")

  for (k in 1:length(Y_list_test)) {

    Y_list_test[[k]] = Y_list_test[[k]][, sort(intersect(names(map), colnames(Y_list_test[[k]]))), drop = FALSE]

  }

  indices = lapply(Y_list_test, function(k) map[colnames(k)])

  attr(indices, "map") = map
  attr(indices, "responses") = names(map)

  return(list(Y_list = Y_list_test, indices_list = indices))

}

subset_observed_data_univariate = function(Y_list, X_list, indices_list, index) {

  new_Y_list = list()
  new_X_list = list()

  for (k in 1:length(indices_list)) {

    if (index %in% indices_list[[k]]) {

      new_Y_list = c(new_Y_list, list(Y_list[[k]][, which(indices_list[[k]] == index), drop = TRUE]))
      new_X_list = c(new_X_list, list(X_list[[k]]))

    }

  }

  Y = do.call(c, new_Y_list)
  X = do.call(rbind, new_X_list)

  return(list(Y = Y, X = X))

}
