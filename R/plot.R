#' Plot actual versus predicted values for each response variable
#'
#' @param Y_list
#' @param X_list
#' @param indices_list
#' @param Beta
#' @param annotations
#'
#' @import ggplot2
#'
#' @return
#' @export
plot_actual_vs_predicted = function(Y_list, X_list, indices_list, Beta, annotations = NULL, order = FALSE) {

  if (ncol(X_list[[1]]) + 1 == nrow(Beta)) {
    X_list = lapply(X_list, function(x) cbind(1, x))
  }

  for (i in 1:length(Y_list)) {
    if (is.null(colnames(Y_list[[i]]))) {
      colnames(Y_list[[i]]) = as.character(indices_list[[i]])
    }
  }

  if (is.null(names(Y_list))) {
    names(Y_list) = paste0("dataset_", 1:length(Y_list))
  }

  if (is.null(colnames(Beta))) {
    colnames(Beta) = as.character(1:ncol(Beta))
  }

  if (!is.matrix(Beta)) {
    Beta = as.matrix(Beta)
  }

  if (!is.null(annotations)) {
    for (i in 1:length(annotations)) {
      if (is.null(names(annotations[[i]]))) names(annotations[[i]]) = colnames(Beta)
    }
    text = apply(do.call(cbind, lapply(names(annotations), function(x) paste0(" ", x, ": ", round(annotations[[x]], 3)))), 1, function(y) paste(y, collapse = "\n"))
    annotation_data = data.frame(text = text, response = names(annotations[[1]]))
  }

  predicted = mapply(X_list, indices_list, FUN = function(x, y) {
    pred = x %*% Beta[, y, drop = FALSE]
    colnames(pred) = colnames(Beta)[y, drop = FALSE]
    pred
  }, SIMPLIFY = FALSE)

  plot_data = do.call(rbind, mapply(Y_list, predicted, names(Y_list), FUN = function(true, pred, name) {
    data = do.call(rbind, lapply(colnames(true), function(x) data.frame(true = true[, x], pred = pred[, x], response = x)))
    data$dataset = name
    data
  }, SIMPLIFY = FALSE))

  if (order & !is.null(annotations)) {
    levels = names(annotations[[1]])[order(annotations[[1]], decreasing = TRUE)]
    plot_data$response = factor(plot_data$response, levels = levels)
    annotation_data$response = factor(annotation_data$response, levels = levels)
  }

  plot = ggplot(plot_data, aes(x = true, y = pred, color = dataset)) + geom_point(size = 0.25) + facet_wrap(~response, scales = "free") + theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside") + labs(x = "actual", y = "predicted")

  if (!is.null(annotations)) {
    plot = plot + geom_text(data = annotation_data, aes(x = -Inf, y = Inf, label = text, vjust = 1, hjust = 0), inherit.aes = FALSE, parse = FALSE, size = 3)
  }

  return(plot)

}

plot_L = function(L_list, metadata_list, color, columns = 1:2) {

  stopifnot(length(color) == 1)
  stopifnot(length(columns) == 2)

  if (is.null(names(L_list))) {
    names(L_list) = paste0("dataset_", 1:length(L_list))
  }

  columns = 3:4

  plot_data = do.call(rbind, mapply(L_list, metadata_list, names(L_list), FUN = function(L, metadata, name) {
    data = data.frame(L$u[, columns] %*% diag(L$d[columns]), metadata[, color], name)
    colnames(data) = c(paste0("L_", columns), color, "dataset")
    data
  }, SIMPLIFY = FALSE))

  ggplot(plot_data, aes_string(x = colnames(plot_data)[1], y = colnames(plot_data)[2], color = colnames(plot_data)[3])) + geom_point(size = 0.25) + facet_wrap(~dataset, scales = "free") + theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside") + facet_wrap(~dataset, ncol = 3, scales = "fixed")


  ggplot(plot_data, aes_string(x = colnames(plot_data)[1], y = colnames(plot_data)[2], color = colnames(plot_data)[3])) + geom_point(size = 0.25) + facet_wrap(~dataset, scales = "free") + theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside")

}

plot_L_density = function(L_list, metadata_list, color, columns = 1:5) {

  stopifnot(length(color) == 1)

  if (is.null(names(L_list))) {
    names(L_list) = paste0("dataset_", 1:length(L_list))
  }

  plot_data = do.call(rbind, lapply(columns, function(column) do.call(rbind, mapply(L_list, metadata_list, names(L_list), FUN = function(L, metadata, name) {
    data = data.frame(L$u[, column], paste0("L_", column), metadata[, color], name)
    colnames(data) = c("L", "column", "color", "dataset")
    data
  }, SIMPLIFY = FALSE))))

  ggplot(plot_data, aes_string(x = "L", y = "column", color = "column")) + ggridges::geom_density_ridges2() + facet_grid(color~dataset) + theme_classic() +
    theme(strip.background = element_blank(), strip.placement = "outside")

}