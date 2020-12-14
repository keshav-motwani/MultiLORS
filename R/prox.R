# l1_prox = function(matrix, lambda, intercept = TRUE) {
#
#   thresholded = matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
#
#   thresholded[matrix > lambda] = matrix[matrix > lambda] - lambda
#
#   thresholded[matrix < -1 * lambda] = matrix[matrix < -1 * lambda] + lambda
#
#   if (intercept) thresholded[1, ] = matrix[1, ]
#
#   return(thresholded)
#
# }
#
# nuclear_prox = function(matrix, gamma) {
#
#   decomp = svd(matrix)
#
#   thresholded = decomp$u %*% tcrossprod(diag(pmax(decomp$d - gamma, 0)), decomp$v)
#
#   return(thresholded)
#
# }
