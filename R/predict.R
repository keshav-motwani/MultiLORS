#' Predict Y
#'
#' @param X
#' @param Beta
#'
#' @return
#' @export
predict_Y = function(X, Beta) {

  return(cbind(1, X) %*% Beta)

}