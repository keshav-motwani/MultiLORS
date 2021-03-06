library(MultiLORS)

result_path = "results/applications/hao_3_prime/"

data = readRDS(file.path(result_path, "data.rds"))
Y_list = data$Y_list
X_list = data$X_list
metadata_list = data$metadata_list

train_data = 1:12
test_data = 13:24

prepared_train = prepare_Y_and_indices_train(Y_list[train_data])
prepared_test = prepare_Y_and_indices_test(Y_list[test_data], prepared_train$indices_list)

MultiLORS_fit = MultiLORS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                          verbose = 1, n_iter = 500, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = FALSE,
                          return_L = FALSE, n_cores = 20, gamma_min_ratio = 1e-4, extra_iter = 5000, extra_iter_threshold = 0.25, standardize_Y = TRUE)

MultiLORS_fit = readRDS(file.path(result_path, "MultiLORS_fit.rds"))

MultiLORS_refit = refit_MultiLORS(MultiLORS_fit, prepared_train$Y_list, X_list[train_data], prepared_train$indices_list, n_cores = 10)

saveRDS(MultiLORS_refit, file.path(result_path, "MultiLORS_refit_fit.rds"))