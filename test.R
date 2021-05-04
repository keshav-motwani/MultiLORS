library(MultiLORS)

set.seed(111)

n_k = rep(1000, 8)
p = 500
q = 30
theta = rep(50, 8)
Sigma = diag(1, nrow = q, ncol = q)
sparsity = 0.9
R2 = 0.9

X_star = simulate_X_list(n_k, p)
L = simulate_L_list(n_k, q, theta)
E = simulate_E_list(n_k, Sigma)
Beta = simulate_Beta(p, q, sparsity, R2, X_star, L, E)
Beta[1, ] = 0

compute_U = function(L, Beta) {

  pinv = svd(tcrossprod(Beta))
  d = 1/pinv$d
  d[pinv$d < 1e-12] = 0

  L %*% t(Beta) %*% (pinv$u %*% diag(d) %*% t(pinv$v))

}

compute_U_list = function(L_list, Beta) {
  lapply(L_list, compute_U, Beta)
}

U = compute_U_list(L, Beta[-1, ])
X = mapply(`+`, X_star, U, SIMPLIFY = FALSE)
Y = compute_Y_list(X_star, Beta, rep(0, length(X_star)), E)
indices = replicate(8, 1:q, simplify = FALSE)

X_star_val = simulate_X_list(n_k, p)
L_val = simulate_L_list(n_k, q, theta)
E_val = simulate_E_list(n_k, Sigma)
U_val = compute_U_list(L_val, Beta[-1, ])
X_val = mapply(`+`, X_star_val, U_val, SIMPLIFY = FALSE)
Y_val = compute_Y_list(X_star_val, Beta, rep(0, length(X_star_val)), E_val)
indices_val = replicate(8, 1:q, simplify = FALSE)

print(system.time({fit = MultiLORS(Y, X_star, indices, Y_val, X_star_val, indices_val, verbose = FALSE, n_iter = 1000, early_stopping = FALSE, return_L = T)}))
(tuning = which_min(fit$tuning$validation$SSE))
L[[1]][1:5, 1:5]
fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]][1:5, 1:5]
mean((L[[1]] - fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]])^2)
mean(as.matrix(Beta - fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta)^2)
max(fit$tuning$validation$avg_R2)
(tuning = which_min(fit$tuning$validation$SSE[1, ]))
mean(as.matrix(Beta - fit$model_fits[[1]][[tuning[1]]]$Beta)^2)

Y = compute_Y_list(X_star, Beta, L, E)
Y_val = compute_Y_list(X_star_val, Beta, L_val, E_val)
print(system.time({fit = MultiLORS(Y, X_star, indices, Y_val, X_star_val, indices_val, verbose = FALSE, n_iter = 1000, early_stopping = FALSE, return_L = T)}))
(tuning = which_min(fit$tuning$validation$SSE))
L[[1]][1:5, 1:5]
fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]][1:5, 1:5]
mean((L[[1]] - fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]])^2)
mean(as.matrix(Beta - fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta)^2)
max(fit$tuning$validation$avg_R2)
(tuning = which_min(fit$tuning$validation$SSE[1, ]))
mean(as.matrix(Beta - fit$model_fits[[1]][[tuning[1]]]$Beta)^2)

Y = compute_Y_list(X_star, Beta, rep(0, length(X_star)), E)
Y_val = compute_Y_list(X_star_val, Beta, rep(0, length(X_star_val)), E_val)
print(system.time({fit = MultiLORS(Y, X, indices, Y_val, X_val, indices_val, verbose = FALSE, n_iter = 1000, early_stopping = FALSE, return_L = T)}))
(tuning = which_min(fit$tuning$validation$SSE))
L[[1]][1:5, 1:5]
fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]][1:5, 1:5]
mean((L[[1]] + fit$model_fits[[tuning[1]]][[tuning[2]]]$L_list[[1]])^2)
mean(as.matrix(Beta - fit$model_fits[[tuning[1]]][[tuning[2]]]$Beta)^2)
max(fit$tuning$validation$avg_R2)
(tuning = which_min(fit$tuning$validation$SSE[1, ]))
mean(as.matrix(Beta - fit$model_fits[[1]][[tuning[1]]]$Beta)^2)

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
