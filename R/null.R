compute_Y_mean = function(Y_list_train, D_list_train) {

  sapply(1:ncol(D_list_train[[1]]), function(x) mean(subset_observed_data_univariate(Y_list_train, NULL, D_list_train, x)$Y))

}