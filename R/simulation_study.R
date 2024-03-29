#' Generate data for simulation study
#'
#' @param R2
#' @param r
#' @param theta
#' @param sigma_2
#' @param sparsity
#' @param p
#' @param q
#' @param O_k_train
#' @param n_k_train
#' @param O_k_validation
#' @param n_k_validation
#' @param O_k_test
#' @param n_k_test
#'
#' @return
#' @export
generate_simulation_data = function(R2,
                                    r,
                                    theta,
                                    sigma_2,
                                    sparsity,
                                    p,
                                    q,
                                    O_k_train,
                                    n_k_train,
                                    O_k_validation,
                                    n_k_validation,
                                    O_k_test,
                                    n_k_test,
                                    replicate) {

  set.seed(replicate)

  theta = rep(theta, r)

  Sigma = matrix(0, nrow = q, ncol = q)
  diag(Sigma) = sigma_2

  n = list(train = n_k_train,
           validation = n_k_validation,
           test = n_k_test)

  O = list(train = O_k_train,
           validation = O_k_validation,
           test = O_k_test)

  X = lapply(n, function(n_k) simulate_X_list(n_k, p))

  L = lapply(n, function(n_k) simulate_L_list(n_k, q, theta))

  E = lapply(n, function(n_k) simulate_E_list(n_k, Sigma))

  Beta = simulate_Beta(p,
                       q,
                       sparsity,
                       R2,
                       unlist(X, recursive = FALSE),
                       unlist(L, recursive = FALSE),
                       unlist(E, recursive = FALSE))

  Y = mapply(
    X_list = X,
    L_list = L,
    E_list = E,
    MoreArgs = list(Beta = Beta),
    compute_Y_list,
    SIMPLIFY = FALSE
  )

  subsetted = mapply(
    Y_list = Y,
    L_list = L,
    E_list = E,
    O_k = O,
    subset_simulated_data,
    SIMPLIFY = FALSE
  )

  indices_full = 1:q

  names = c("train", "validation", "test")
  output = lapply(names,
                  function(name)
                    list(
                      full = format_simulation_output_full(Y[[name]], X[[name]], L[[name]], E[[name]], indices_full),
                      subsetted = format_simulation_output_subsetted(subsetted[[name]], X[[name]])
                    ))
  names(output) = names

  output$Beta = Beta

  return(output)

}

format_simulation_output_full = function(Y_list, X_list, L_list, E_list, indices) {

  return(list(
    Y_list = Y_list,
    X_list = X_list,
    L_list = L_list,
    E_list = E_list,
    indices_list = replicate(length(Y_list), indices, simplify = FALSE)
  ))

}

format_simulation_output_subsetted = function(subsetted, X_list) {

  return(c(list(X_list = X_list),
           subsetted))

}

#' Expand parameters one at a time given defaults and considered values
#'
#' @param considered_values
#' @param defaults
#' @param n_replicates
#' @param methods
#'
#' @return
#' @export
expand_parameters = function(considered_values,
                             defaults,
                             n_replicates,
                             methods) {

  parameter_list = list()

  for (name in names(considered_values)) {
    values = considered_values[[name]]
    for (replicate in 1:n_replicates) {
      for (value in values) {
        for (method in methods) {
          params = c(
            defaults,
            experiment = name,
            replicate = replicate,
            method = method
          )
          params[name] = value
          parameter_list = c(parameter_list, list(params))
        }
      }
    }
  }

  return(parameter_list)
}

#' Evaluate parameters
#'
#' @param parameters
#'
#' @return
#' @export
evaluate_parameters = function(parameters) {

  method_function = get(paste0("fit_model_", parameters$method))

  data = do.call(generate_simulation_data, parameters[formalArgs(generate_simulation_data)])
  fit = method_function(data)

  performance = compute_performance(data$Beta,
                                    fit$Beta,
                                    data$test$full$Y_list,
                                    data$test$full$X_list,
                                    data$test$full$indices_list,
                                    fit$Y_list_train,
                                    fit$indices_list_train)

  return(list(parameters = parameters, result = performance))

}

fit_model_MultiLORS = function(data) {

  Y_list_train = data$train$subsetted$Y_list
  indices_list_train = data$train$subsetted$indices_list

  fit = MultiLORS(
    Y_list_train,
    data$train$subsetted$X_list,
    indices_list_train,
    data$validation$subsetted$Y_list,
    data$validation$subsetted$X_list,
    data$validation$subsetted$indices_list,
    verbose = 1,
    n_iter = 1000,
    tolerance = 1e-6
  )

  return(list(Beta = fit$best_Beta, Y_list_train= Y_list_train, indices_list_train = indices_list_train))
}

fit_model_glmnet = function(data) {

  Y_list_train = data$train$subsetted$Y_list
  indices_list_train = data$train$subsetted$indices_list

  Beta = fit_glmnet(Y_list_train,
             data$train$subsetted$X_list,
             indices_list_train,
             data$validation$subsetted$Y_list,
             data$validation$subsetted$X_list,
             data$validation$subsetted$indices_list)$best_Beta

  return(list(Beta = Beta, Y_list_train= Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_ALL_MultiLORS = function(data) {

  Y_list_train = data$train$full$Y_list
  indices_list_train = data$train$full$indices_list

  Beta = MultiLORS(
    Y_list_train,
    data$train$full$X_list,
    indices_list_train,
    data$validation$full$Y_list,
    data$validation$full$X_list,
    data$validation$full$indices_list,
    verbose = 1,
    n_iter = 1000,
    tolerance = 1e-6
  )$best_Beta

  return(list(Beta = Beta, Y_list_train= Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_L_glmnet = function(data) {

  Y_list_train = mapply(x = data$train$subsetted$Y_list, y = data$train$subsetted$L_list, function(x, y) x - y, SIMPLIFY = FALSE)
  indices_list_train = data$train$subsetted$indices_list

  Beta = fit_glmnet(Y_list_train,
                    data$train$subsetted$X_list,
                    indices_list_train,
             mapply(x = data$validation$subsetted$Y_list, y = data$validation$subsetted$L_list, function(x, y) x - y, SIMPLIFY = FALSE),
             data$validation$subsetted$X_list,
             data$validation$subsetted$indices_list)$best_Beta

  return(list(Beta = Beta, Y_list_train= Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_L_MultiLORS = function(data) {

  Y_list_train = mapply(x = data$train$subsetted$Y_list, y = data$train$subsetted$L_list, function(x, y) x - y, SIMPLIFY = FALSE)
  indices_list_train = data$train$subsetted$indices_list

  Beta = MultiLORS(
    Y_list_train,
    data$train$subsetted$X_list,
    indices_list_train,
    mapply(x = data$validation$subsetted$Y_list, y = data$validation$subsetted$L_list, function(x, y) x - y, SIMPLIFY = FALSE),
    data$validation$subsetted$X_list,
    data$validation$subsetted$indices_list,
    verbose = 1,
    n_iter = 1000,
    tolerance = 1e-6
  )$best_Beta

  return(list(Beta = Beta, Y_list_train = Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_ALL_glmnet = function(data) {

  Y_list_train = data$train$full$Y_list
  indices_list_train = data$train$full$indices_list

  Beta = fit_glmnet(Y_list_train,
                    data$train$full$X_list,
                    indices_list_train,
                    data$validation$full$Y_list,
                    data$validation$full$X_list,
                    data$validation$full$indices_list)$best_Beta

  return(list(Beta = Beta, Y_list_train = Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_L_ALL_glmnet = function(data) {

  Y_list_train = mapply(x = data$train$full$Y_list, y = data$train$full$L_list, function(x, y) x - y, SIMPLIFY = FALSE)
  indices_list_train = data$train$full$indices_list

  Beta = fit_glmnet(Y_list_train,
             data$train$full$X_list,
             indices_list_train,
             mapply(x = data$validation$full$Y_list, y = data$validation$full$L_list, function(x, y) x - y, SIMPLIFY = FALSE),
             data$validation$full$X_list,
             data$validation$full$indices_list)$best_Beta

  return(list(Beta = Beta, Y_list_train = Y_list_train, indices_list_train = indices_list_train))

}

fit_model_ORC_L_ALL_MultiLORS = function(data) {

  Y_list_train = mapply(x = data$train$full$Y_list, y = data$train$full$L_list, function(x, y) x - y, SIMPLIFY = FALSE)
  indices_list_train = data$train$full$indices_list

  Beta = MultiLORS(
    Y_list_train,
    data$train$full$X_list,
    indices_list_train,
    mapply(x = data$validation$full$Y_list, y = data$validation$full$L_list, function(x, y) x - y, SIMPLIFY = FALSE),
    data$validation$full$X_list,
    data$validation$full$indices_list,
    verbose = 1,
    n_iter = 1000,
    tolerance = 1e-6
  )$best_Beta

  return(list(Beta = Beta, Y_list_train = Y_list_train, indices_list_train = indices_list_train))

}

compute_performance = function(Beta, Beta_hat, Y_list_test, X_list_test, indices_list_test, Y_list_train, indices_list_train) {

  Beta_SSE = sum((Beta - Beta_hat)^2)

  SST = sum(compute_SST(Y_list_test, indices_list_test, Y_list_train, indices_list_train))
  SSE = sum(compute_SSE(Y_list_test, X_list_test, indices_list_test, Beta_hat))
  R2 = 1 - SSE/SST

  return(list(Beta_SSE = Beta_SSE, test_R2 = R2, test_SSE = SSE, test_SST = SST, Beta_hat = Beta_hat, Beta = Beta))

}
