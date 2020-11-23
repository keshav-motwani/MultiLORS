fit_glmnet = function(Y_list,
                      X_list,
                      D_list,
                      Y_list_validation,
                      X_list_validation,
                      D_list_validation) {

  Beta = matrix(ncol = ncol(D_list[[1]]), nrow = ncol(X_list[[1]]) + 1)

  for (index in 1:ncol(D_list[[1]])) {

    data = subset_observed_data_univariate(Y_list, X_list, D_list, index)
    data_validation = subset_observed_data_univariate(Y_list_validation, X_list_validation, D_list_validation, index)

    fit = fit_glmnet_univariate(data$Y, data$X, data_validation$Y, data_validation$X)
    Beta[, index] = fit

  }

  return(Beta)

}

fit_glmnet_univariate = function(Y, X, Y_validation, X_validation) {

  fit = glmnet::glmnet(X, Y)

  errors = numeric(length = length(fit$lambda))

  predictions = predict(fit, newx = X_validation)

  for (i in 1:length(fit$lambda)) {

    errors[i] = sum((Y_validation - predictions[, i]) ^ 2)

  }

  Beta = coef(fit)[, which.min(errors)]

  return(Beta)

}
