library(MultiLORS)

set.seed(10, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

n_k = rep(200, 5)
p = 50
q = 5
theta = rep(100, 2)
Sigma = diag(1, nrow = q, ncol = q)
sparsity = 0.9
R2 = 0.75
mean = 5

X = simulate_X_list(n_k, p, mean)
L = simulate_L_list(n_k, q, theta)
E = simulate_E_list(n_k, Sigma)
Beta = simulate_Beta(p, q, sparsity, R2, X, L, E)
print(Beta[1:5, 1:5])
Y = compute_Y_list(X, Beta, L, E)
indices = replicate(5, 1:q, simplify = FALSE)

X_val = simulate_X_list(n_k, p, mean)
L_val = simulate_L_list(n_k, q, theta)
E_val = simulate_E_list(n_k, Sigma)
Y_val = compute_Y_list(X_val, Beta, L_val, E_val)
indices_val = replicate(5, 1:q, simplify = FALSE)

X_test = simulate_X_list(n_k, p, mean)
L_test = simulate_L_list(n_k, q, theta)
E_test = simulate_E_list(n_k, Sigma)
Y_test = compute_Y_list(X_test, Beta, L_test, E_test)
indices_test = replicate(5, 1:q, simplify = FALSE)

set.seed(Sys.time()) # results should be same regardless of seed

print(system.time({MultiLORS_fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = 0, n_iter = 250, n_cores = 4, early_stopping = FALSE)}))
MultiLORS_refit = refit_MultiLORS(MultiLORS_fit, Y, X, indices, Y_val, X_val, indices_val)
glmnet_fit = fit_glmnet(Y, X, indices, Y_val, X_val, indices_val)
glmnet_refit = refit_glmnet(glmnet_fit, Y, X, indices, Y_val, X_val, indices_val)
OLS_fit = fit_OLS(Y, X, indices)

best_indices_MultiLORS = which_min(MultiLORS_fit$tuning$validation$SSE)
best_indices_MultiLORS_refit = which_min(MultiLORS_refit$tuning$validation$SSE)
best_index_glmnet = which_min(glmnet_fit$tuning$validation$SSE)
best_index_glmnet_refit = which_min(glmnet_refit$tuning$validation$SSE)

print(best_indices_MultiLORS)
print(best_indices_MultiLORS_refit)
print(best_index_glmnet)
print(best_index_glmnet_refit)

MultiLORS_Beta_hat = MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta
MultiLORS_refit_Beta_hat = MultiLORS_refit$model_fits[[best_indices_MultiLORS_refit[1]]][[best_indices_MultiLORS_refit[2]]]$Beta
glmnet_Beta_hat = glmnet_fit$model_fits[[best_index_glmnet]]$Beta
glmnet_refit_Beta_hat = glmnet_refit$model_fits[[best_index_glmnet_refit]]$Beta
OLS_Beta_hat = OLS_fit

sum((Beta[-1, ] - as.matrix(MultiLORS_Beta_hat)[-1, ]) ^ 2)
sum((Beta[-1, ] - as.matrix(MultiLORS_refit_Beta_hat)[-1, ]) ^ 2)
sum((Beta[-1, ] - as.matrix(glmnet_Beta_hat)[-1, ]) ^ 2)
sum((Beta[-1, ] - as.matrix(glmnet_refit_Beta_hat)[-1, ]) ^ 2)
sum((Beta[-1, ] - as.matrix(OLS_Beta_hat)[-1, ]) ^ 2)

compute_avg_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat)
compute_avg_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_refit_Beta_hat)
compute_avg_R2(Y_test, X_test, indices_test, Y, indices, glmnet_Beta_hat)
compute_avg_R2(Y_test, X_test, indices_test, Y, indices, glmnet_refit_Beta_hat)
compute_avg_R2(Y_test, X_test, indices_test, Y, indices, OLS_Beta_hat)

#
# compute_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat)
# compute_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat)
# compute_SSE(Y_test, X_test, indices_test, MultiLORS_Beta_hat)
#

# compute_avg_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat)
# compute_avg_correlation(Y_test, X_test, indices_test, glmnet_Beta_hat)
#
# MultiLORS_performance = list("R2" = compute_R2(Y_test, X_test, indices_test, Y, indices, MultiLORS_Beta_hat),
#                              "rho" = compute_correlation(Y_test, X_test, indices_test, MultiLORS_Beta_hat))
# plot_actual_vs_predicted(Y_test, X_test, indices_test, MultiLORS_Beta_hat, MultiLORS_performance)
#
# glmnet_performance = list("R2" = compute_R2(Y_test, X_test, indices_test, Y, indices, glmnet_Beta_hat),
#                              "rho" = compute_correlation(Y_test, X_test, indices_test, glmnet_Beta_hat))
# plot_actual_vs_predicted(Y_test, X_test, indices_test, glmnet_Beta_hat, glmnet_performance)
