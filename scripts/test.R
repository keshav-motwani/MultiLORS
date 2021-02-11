library(MultiLORS)

set.seed(1)

n_k = rep(500, 3)
p = 50
q = 10
theta = rep(100, 3)
Sigma = diag(1, nrow = q, ncol = q)
sparsity = 0.9
R2 = 0.9
mean = 5

X = simulate_X_list(n_k, p, mean)
L = simulate_L_list(n_k, q, theta)
E = simulate_E_list(n_k, Sigma)
Beta = simulate_Beta(p, q, sparsity, R2, X, L, E)
Y = compute_Y_list(X, Beta, L, E)
indices = list(1:q, 1:q, 1:q)

X_val = simulate_X_list(n_k, p, mean)
L_val = simulate_L_list(n_k, q, theta)
E_val = simulate_E_list(n_k, Sigma)
Y_val = compute_Y_list(X_val, Beta, L_val, E_val)
indices_val = list(1:q, 1:q, 1:q)

X_test = simulate_X_list(n_k, p, mean)
L_test = simulate_L_list(n_k, q, theta)
E_test = simulate_E_list(n_k, Sigma)
Y_test = compute_Y_list(X_test, Beta, L_test, E_test)
indices_test = list(1:q, 1:q, 1:q)

print(system.time({MultiLORS_fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = 0, n_iter = 100, n_cores = 4, early_stopping = FALSE)}))
print(MultiLORS_fit$tuning$validation$n_iter)
glmnet_fit = fit_glmnet(Y, X, indices, Y_val, X_val, indices_val)

best_indices_MultiLORS = which_min(MultiLORS_fit$tuning$validation$SSE)
best_index_glmnet = which_min(glmnet_fit$tuning$validation$SSE)
sum((Beta - as.matrix(MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta)) ^ 2)
sum((Beta - as.matrix(glmnet_fit$model_fits[[best_index_glmnet]]$Beta)) ^ 2)

MultiLORS_Beta_hat = MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta
glmnet_Beta_hat = glmnet_fit$model_fits[[best_index_glmnet]]$Beta

compute_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat)
compute_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat)
compute_SSE(Y_test, X_test, indices_test, MultiLORS_Beta_hat)

compute_avg_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat)
compute_avg_R2(Y_test, X_test, indices_test, Y, indices, glmnet_Beta_hat)
compute_avg_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat)
compute_avg_correlation(Y_test, X_test, indices_test, glmnet_Beta_hat)

MultiLORS_performance = list("R2" = compute_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat),
                             "rho" = compute_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat))
plot_actual_vs_predicted(Y_test, X_test, indices_test, MultiLORS_Beta_hat, MultiLORS_performance)

glmnet_performance = list("R2" = compute_R2(Y_test, X_test, indices_test, Y, indices, glmnet_Beta_hat),
                             "rho" = compute_correlation(Y_test, X_test, indices_test, glmnet_Beta_hat))
plot_actual_vs_predicted(Y_test, X_test, indices_test, glmnet_Beta_hat, glmnet_performance)
