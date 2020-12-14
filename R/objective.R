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

error = function(Y, Y_hat) {

  sum((Y - Y_hat) ^ 2)

}

sst = function(Y_train_means, Y_list_test, D_list_test) {

  sst = sum(sapply(1:length(Y_list_test), function(k) error(Y_list_test[[k]], tcrossprod(rep(1, nrow(Y_list_test[[k]])), Y_train_means[diag(D_list_test[[k]]) == 1]))))

  return(sst)

}
