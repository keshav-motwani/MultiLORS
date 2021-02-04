# evaluate_objective = function(Y_list, X_list, L_list, indices_list, Beta, lambda, gamma, gamma_weights) {
#
#   p1 = evaluate_g(Y_list, X_list, L_list, indices_list, Beta)
#
#   p2 = nuclear_norm_penalty(L_list, gamma, gamma_weights)
#
#   p3 = l1_penalty(Beta, lambda)
#
#   return(p1 + p2 + p3)
#
# }
#
# evaluate_g = function(Y_list, X_list, L_list, indices_list, Beta)  {
#
#   p1 = 0
#
#   for (k in 1:length(Y_list)) {
#
#     observed_indices = indices_list[[k]]
#
#     p1 = p1 + error(Y_list[[k]], (X_list[[k]] %*% Beta[, observed_indices]) + L_list[[k]])
#
#   }
#
#   return(p1 / 2)
#
# }
#
# l1_penalty = function(Beta, lambda) {
#
#   return(lambda * sum(abs(Beta[-1, ])))
#
# }
#
# nuclear_norm_penalty = function(L_list, gamma, gamma_weights) {
#
#   value = 0
#
#   for (k in 1:length(L_list)) {
#
#     gamma_k = gamma_weights[[k]] * gamma
#
#     decomp = svd(L_list[[k]])
#
#     value = value + sum(decomp$d) * gamma_k
#
#   }
#
#   return(value)
#
# }
#
# compute_error = function(Y_list, X_list, indices_list, Beta) {
#
#   error = 0
#
#   for (k in 1:length(Y_list)) {
#
#     error = error + sum((Y_list[[k]] - X_list[[k]] %*% Beta[, indices_list[[k]]]) ^ 2)
#
#   }
#
#   return(error)
#
# }
