#' simulate E error matrix with specified covariance matrix
#'
#' @param n
#' @param Sigma
#'
#' @return
#' @export
simulate_E = function(n, Sigma) {

  eo = eigen(Sigma)

  E = matrix(rnorm(n * nrow(Sigma)), nrow = n) %*% eo$vec %*% diag(eo$val^.5) %*% t(eo$vec)

  return(E)

}

#' simulate list of E error matrices with specified covariance matrix
#'
#' @param n_k
#' @param Sigma
#'
#' @return
#' @export
simulate_E_list = function(n_k, Sigma) {

  lapply(n_k, function(n) simulate_E(n, Sigma))

}

#' simulate L batch effect matrix
#'
#' @param n
#' @param q
#' @param r
#' @param variance
#'
#' @return
#' @export
simulate_L = function(n, q, theta) {

  r = length(theta)

  R = matrix(rnorm(q * q), ncol = q)

  Q = svd(R)$v[, 1:r]

  Sigma_sqrt = tcrossprod(Q %*% diag(sqrt(theta)), Q)

  W = matrix(rnorm(n * q), ncol = q)

  L = W %*% Sigma_sqrt

  return(L)

}

#' simulate list of L batch effect matrices
#'
#' @param n_k
#' @param q
#' @param r
#' @param variance
#'
#' @return
#' @export
simulate_L_list = function(n_k, q, theta) {

  lapply(n_k, function(n) simulate_L(n, q, theta))

}

solve_c = function(X, Beta, L, E, R2, i) {
  min_function = function(c) {
    abs(((var(c * X %*% Beta[-1, i] + E[, i]) - var(-L[, i] + E[, i])) / var(-L[, i] + E[, i])) - (R2 / (1 - R2)))
  }
  result = optimize(min_function, c(0, 20), tol = 1e-100)
  print(result)
  result$minimum
}

#' simulate row-sparse Beta matrix
#'
#' @param p
#' @param q
#'
#' @return
#' @export
simulate_Beta = function(p, q, sparsity, R2 = NULL, X = NULL, L = NULL, E = NULL) {

  if (is.list(X)) X = do.call(rbind, X)
  if (is.list(L)) L = do.call(rbind, L)
  if (is.list(E)) E = do.call(rbind, E)

  Beta = matrix(runif(p * q, min = -3, max = 3), nrow = p, ncol = q)

  Beta[sample(1:(p*q), size = p * q * sparsity)] = 0

  Beta = rbind(runif(q, min = -10, max = 10), Beta)

  if (!is.null(R2)) {

    Beta[-1, ] = Beta[-1, ] %*% diag(1/sqrt(colSums(Beta[-1, ] ^ 2)))

    c = sapply(1:q, function(i) solve_c(X, Beta, L, E, R2, i))

    Beta[-1, ] = Beta[-1, ] %*% diag(c)

  }

  Y = cbind(1, X) %*% Beta + E
  SST = colSums((Y - (rep(1, nrow(Y)) %*% crossprod(rep(1/nrow(Y), nrow(Y)), Y))) ^ 2)
  SSE = colSums((Y - cbind(1, X) %*% Beta - L) ^ 2)
  print((SST - SSE) / SST)

  return(Beta)

}

#' simulate X matrix
#'
#' @param n
#' @param p
#'
#' @return
#' @export
simulate_X = function(n, p, mean = 0) {

  Sigma = matrix(nrow = p, ncol = p)

  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] = 0.9 ^ abs(i - j)
    }
  }

  eo = eigen(Sigma)

  X = matrix(rnorm(n * nrow(Sigma)), nrow = n) %*% eo$vec %*% diag(eo$val^.5) %*% t(eo$vec)

  X = X + mean

  return(X)

}

#' simulate list of X matrices
#'
#' @param n_k
#' @param p
#'
#' @return
#' @export
simulate_X_list = function(n_k, p, mean = 0) {

  return(lapply(n_k, function(n) simulate_X(n, p, mean)))

}

#' Compute Y given X, Beta, L, E
#'
#' @param X
#' @param Beta
#' @param L
#' @param E
#'
#' @return
#' @export
compute_Y = function(X, Beta, L, E) {

  Y = predict_Y(X, Beta) + L + E

  return(Y)

}

#' Compute list of Y matrices given X, Beta, L, E
#'
#' @param X_list
#' @param Beta
#' @param L_list
#' @param E_list
#'
#' @return
#' @export
compute_Y_list = function(X_list, Beta, L_list, E_list) {

  mapply(compute_Y, X = X_list, L = L_list, E = E_list, MoreArgs = list(Beta = Beta), SIMPLIFY = FALSE)

}

#' Subset simulated data to have different observed responses in each dataset
#'
#' @param Y_list
#' @param L_list
#' @param E_list
#' @param O_k
#'
#' @return
#' @export
subset_simulated_data = function(Y_list, L_list, E_list, O_k) {

  if (!is.list(O_k)) {

    q = ncol(Y_list[[1]])

    total_remaining = sum(O_k)
    available_indices = 1:q

    indices = vector("list", length = length(O_k))

    while(total_remaining > 0) {

      for (k in 1:length(O_k)) {

        if (length(indices[[k]]) < O_k[k]) {

          result = sample_state(available_indices, indices[[k]], 1:q)
          indices[[k]] = sort(c(indices[[k]], result$sampled))
          available_indices = result$values

          total_remaining = total_remaining - 1

        }

      }

    }

  } else {

    indices = O_k

  }

  indices = lapply(indices, function(x) sort(unique(x)))

  Y_list = mapply(function(m, i) m[, i, drop = FALSE], m = Y_list, i = indices, SIMPLIFY = FALSE)
  L_list = mapply(function(m, i) m[, i, drop = FALSE], m = L_list, i = indices, SIMPLIFY = FALSE)
  E_list = mapply(function(m, i) m[, i, drop = FALSE], m = E_list, i = indices, SIMPLIFY = FALSE)

  return(list(Y_list = Y_list,
         L_list = L_list,
         E_list = E_list,
         indices_list = indices))

}

sample_state = function(available_values, used_values, all_values) {

  if (length(available_values) == 0) available_values = all_values

  if (length(setdiff(available_values, used_values)) == 0) {
    sampled = sample_1(setdiff(all_values, used_values))
  } else {
    sampled = sample_1(setdiff(available_values, used_values))
  }

  available_values = available_values[available_values != sampled]

  return(list(sampled = sampled, values = available_values))

}

sample_1 = function(x) {

  if (length(x) == 1) return(x)
  else return(sample(x, 1))

}