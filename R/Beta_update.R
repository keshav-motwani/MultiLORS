# update_Beta = function(Y_list, X_list, L_list, q, indices_list, XtX_list, XtY_list, Beta_old, lambda, s_Beta, s) {
#
#   gradient = compute_gradient_Beta(X_list = X_list,
#                                    L_list = L_list,
#                                    q = q,
#                                    indices_list = indices_list,
#                                    XtX_list = XtX_list,
#                                    XtY_list = XtY_list,
#                                    Beta_old = Beta_old)
#
#   line_search = TRUE
#   shrinkage = 0.5
#
#   # https://people.eecs.berkeley.edu/~elghaoui/Teaching/EE227A/lecture18.pdf
#
#     g_old = evaluate_g(Y_list, X_list, L_list, indices_list, Beta_old)
#
#     while(line_search & s != s_Beta) {
#
#       Beta = l1_prox(Beta_old - (s * gradient), s * lambda)
#
#       g_new = evaluate_g(Y_list, X_list, L_list, indices_list, Beta)
#
#       difference = (Beta_old - Beta) / s
#
#       if (g_new > g_old - s * sum(gradient * difference) + 0.5 * s * sum(difference^2)) {
#
#         s = max(shrinkage * s, s_Beta)
#
#       } else {
#
#         line_search = FALSE
#
#       }
#
#     }
#
#   # print(paste0("s: ", s, "; s_Beta: ", s_Beta))
#
#   return(Beta)
#
# }
#
# compute_gradient_Beta = function(X_list, L_list, q, indices_list, XtX_list, XtY_list, Beta_old) {
#
#   gradient = matrix(0, nrow = ncol(X_list[[1]]), ncol = q)
#
#   for (k in 1:length(X_list)) {
#
#     one = XtX_list[[k]] %*% Beta_old[, indices_list[[k]], drop = FALSE]
#
#     two = crossprod(X_list[[k]], L_list[[k]])
#
#     three = XtY_list[[k]]
#
#     gradient[, indices_list[[k]]] = gradient[, indices_list[[k]]] + one + two - three
#
#   }
#
#   return(gradient)
#
# }
#
# compute_s_Beta = function(XtX_list, p, q, dataset_indices_list) {
#
#   max_eigenvalues = numeric(q)
#
#   for (j in 1:length(max_eigenvalues)) {
#
#     P_j = dataset_indices_list[[j]]
#
#     block = Reduce('+', XtX_list[P_j])
#
#     eigenvalues = eigen(block, TRUE)$values
#
#     max_eigenvalues[j] = eigenvalues[1]
#
#   }
#
#   return(1 / max(max_eigenvalues))
#
# }
