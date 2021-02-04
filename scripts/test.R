library(MultiLORS)

set.seed(1)
future::plan("multisession", workers = 4)

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

print(system.time({MultiLORS_fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = 0, n_iter = 100, n_cores = 4)}))

glmnet_fit = fit_glmnet(Y, X, indices, Y_val, X_val, indices_val)

best_indices_MultiLORS = which_min(MultiLORS_fit$validation_error)
best_index_glmnet = which_min(glmnet_fit$validation_error)
sum((Beta - as.matrix(MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta)) ^ 2)
sum((Beta - as.matrix(glmnet_fit$model_fits[[best_index_glmnet]]$Beta)) ^ 2)

Beta_hat = MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta

compute_correlation(Y_test, X_test, indices_test, Beta_hat)
compute_R2(Y_test, X_test, indices_test, Y, indices, Beta_hat)
compute_SSE(Y_test, X_test, indices_test, Beta_hat)

compute_avg_R2(Y_test, X_test, indices_test, Y, indices, Beta_hat)
compute_weighted_avg_R2(Y_test, X_test, indices_test, Y, indices, Beta_hat)
compute_avg_correlation(Y_test, X_test, indices_test, Beta_hat)

annotations = list("rho" = compute_correlation(Y_test, X_test, indices_test, Beta_hat),
                   "R^2" = compute_R2(Y_test, X_test, indices_test, Y, indices, Beta_hat))
plot_actual_vs_predicted(Y_test, X_test, indices_test, Beta_hat, annotations)
