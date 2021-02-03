library(MultiLORS)

set.seed(1)

n_k = rep(500, 3)
p = 50
q = 10
theta = rep(100, 3)
Sigma = diag(1, nrow = q, ncol = q)
sparsity = 0.9
R2 = 0.9

X = simulate_X_list(n_k, p)
L = simulate_L_list(n_k, q, theta)
E = simulate_E_list(n_k, Sigma)
Beta = simulate_Beta(p, q, sparsity, R2, X, L, E)
Y = compute_Y_list(X, Beta, L, E)
indices = list(1:q, 1:q, 1:q)

X_val = simulate_X_list(n_k, p)
L_val = simulate_L_list(n_k, q, theta)
E_val = simulate_E_list(n_k, Sigma)
Y_val = compute_Y_list(X_val, Beta, L_val, E_val)
indices_val = list(1:q, 1:q, 1:q)

X_test = simulate_X_list(n_k, p)
L_test = simulate_L_list(n_k, q, theta)
E_test = simulate_E_list(n_k, Sigma)
Y_test = compute_Y_list(X_test, Beta, L_test, E_test)
indices_test = list(1:q, 1:q, 1:q)

print(system.time({fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = FALSE, n_iter = 100, n_cores = 4)}))

print(system.time({fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = FALSE, n_iter = 100)}))


# tuning = which(fit$weighted_avg_validation_R2 == max(fit$weighted_avg_validation_R2, na.rm = TRUE), arr.ind = TRUE)

# sum((as.matrix(fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta) - Beta) ^ 2)

# compute_SSE(Y_test, X_test, indices_test, as.matrix(fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta))
# compute_R2(Y_test, X_test, indices_test, Y, indices, as.matrix(fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta))

# fit_L0 = fit_glmnet(Y, X, indices, Y_val, X_val, indices_val)

# tuning = which.max(fit_L0$weighted_avg_validation_R2)

# sum((as.matrix(fit_L0$model_fits[[tuning]]$Beta) - Beta) ^ 2)

# subset = subset_simulated_data(Y, L, E, c(4, 4, 4))
# Y = subset$Y_list
# L = subset$L_list
# E = subset$E_list
# indices = subset$indices_list

# subset_val = subset_simulated_data(Y_val, L_val, E_val, c(2, 2, 2))
# Y_val = subset_val$Y_list
# L_val = subset_val$L_list
# E_val = subset_val$E_list
# indices_val = subset_val$indices_list

# system.time({fit_subsetted = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = FALSE, n_iter = 100, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = TRUE)})

# tuning = which(fit_subsetted$weighted_avg_validation_R2 == max(fit_subsetted$weighted_avg_validation_R2, na.rm = TRUE), arr.ind = TRUE)

# sum((as.matrix(fit_subsetted$model_fits[[tuning[1]]][[tuning[2]]]$Beta) - Beta) ^ 2)

# compute_SSE(Y_test, X_test, indices_test, as.matrix(fit_subsetted$model_fits[[tuning[1]]][[tuning[2]]]$Beta))
# compute_R2(Y_test, X_test, indices_test, Y, indices, as.matrix(fit_subsetted$model_fits[[tuning[1]]][[tuning[2]]]$Beta))

# plot(compute_R2(Y_test, X_test, indices_test, Y, indices, as.matrix(fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta)) -
#      compute_R2(Y_test, X_test, indices_test, Y, indices, as.matrix(fit_subsetted$model_fits[[tuning[1]]][[tuning[2]]]$Beta)))

# fit_L0 = fit_glmnet(Y, X, indices, Y_val, X_val, indices_val)

# tuning = which.max(fit_L0$weighted_avg_validation_R2)

# sum((as.matrix(fit_L0$model_fits[[tuning]]$Beta) - Beta) ^ 2)

# compute_R2(Y_test, X_test, indices_test, Y, indices, as.matrix(fit_L0$model_fits[[tuning]]$Beta))
