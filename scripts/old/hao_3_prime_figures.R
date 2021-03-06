library(MultiLORS)
library(ggplot2)

result_path = "results/applications/hao_3_prime/"

load(file.path(result_path, "data.RData"))

MultiLORS_fit = readRDS(file.path(result_path, "MultiLORS_fit.rds"))
glmnet_fit = readRDS(file.path(result_path, "glmnet_fit.rds"))

best_indices_MultiLORS = which_min(MultiLORS_fit$tuning$validation$SSE)
best_index_glmnet = which_min(glmnet_fit$tuning$validation$SSE)

MultiLORS_Beta_hat = MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta
glmnet_Beta_hat = glmnet_fit$model_fits[[best_index_glmnet]]$Beta

compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_Beta_hat)
compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_Beta_hat)

plot_performance_comparison = function(performance_1, performance_2, name_1, name_2, text_size = 3, return_plots = FALSE) {

  plots = list()

  for (metric in names(performance_1)) {

    data = data.frame(performance_1[[metric]], performance_2[[metric]], names(performance_1[[1]]))
    colnames(data) = c(name_1, name_2, "response")

    plot = ggplot(data, aes_string(x = name_1, y = name_2)) + geom_point(size = 0.05, color = "red") + geom_text(aes(label = response), size = text_size) + labs(title = metric) + theme_classic() + geom_abline(slope=1, linetype = "dashed")
    plots = c(plots, list(plot))

  }

  if (return_plots) return(plots)

  patchwork::wrap_plots(plots)

}

MultiLORS_annotations = list("R2" = compute_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_Beta_hat),
                             "rho" = compute_correlation(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, MultiLORS_Beta_hat))

glmnet_annotations = list("R2" = compute_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_Beta_hat),
                          "rho" = compute_correlation(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, glmnet_Beta_hat))

plot_performance_comparison(glmnet_annotations, MultiLORS_annotations, "glmnet", "MultiLORS")
ggsave(file.path(result_path, "comparison.jpeg"), height = 4, width = 8, limitsize = FALSE)

MultiLORS_plot = plot_actual_vs_predicted(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations)
ggsave(file.path(result_path, "MultiLORS_test_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
MultiLORS_plot = plot_actual_vs_predicted(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations, order = TRUE)
ggsave(file.path(result_path, "MultiLORS_test_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)

glmnet_plot = plot_actual_vs_predicted(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, glmnet_Beta_hat, glmnet_annotations)
ggsave(file.path(result_path, "glmnet_test_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
glmnet_plot = plot_actual_vs_predicted(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, glmnet_Beta_hat, glmnet_annotations, order = TRUE)
ggsave(file.path(result_path, "glmnet_test_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)


MultiLORS_annotations = list("R2" = compute_R2(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_Beta_hat),
                             "rho" = compute_correlation(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, MultiLORS_Beta_hat))
MultiLORS_plot = plot_actual_vs_predicted(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations)
ggsave(file.path(result_path, "MultiLORS_train_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
MultiLORS_plot = plot_actual_vs_predicted(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations, order = TRUE)
ggsave(file.path(result_path, "MultiLORS_train_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)

glmnet_annotations = list("R2" = compute_R2(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_Beta_hat),
                          "rho" = compute_correlation(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, glmnet_Beta_hat))
glmnet_plot = plot_actual_vs_predicted(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, glmnet_Beta_hat, glmnet_annotations)
ggsave(file.path(result_path, "glmnet_train_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
glmnet_plot = plot_actual_vs_predicted(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, glmnet_Beta_hat, glmnet_annotations, order = TRUE)
ggsave(file.path(result_path, "glmnet_train_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)


MultiLORS_annotations = list("R2" = compute_R2(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_Beta_hat),
                             "rho" = compute_correlation(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, MultiLORS_Beta_hat))
MultiLORS_plot = plot_actual_vs_predicted(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations)
ggsave(file.path(result_path, "MultiLORS_validation_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
MultiLORS_plot = plot_actual_vs_predicted(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, MultiLORS_Beta_hat, MultiLORS_annotations, order = TRUE)
ggsave(file.path(result_path, "MultiLORS_validation_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)

glmnet_annotations = list("R2" = compute_R2(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_Beta_hat),
                          "rho" = compute_correlation(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, glmnet_Beta_hat))
glmnet_plot = plot_actual_vs_predicted(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, glmnet_Beta_hat, glmnet_annotations)
ggsave(file.path(result_path, "glmnet_validation_set_plot.jpeg"), height = 50, width = 50, limitsize = FALSE)
glmnet_plot = plot_actual_vs_predicted(prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list, glmnet_Beta_hat, glmnet_annotations, order = TRUE)
ggsave(file.path(result_path, "glmnet_validation_set_plot_ordered.jpeg"), height = 50, width = 50, limitsize = FALSE)


