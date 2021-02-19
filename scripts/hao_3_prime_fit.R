library(MultiLORS)
library(CITEseqData)
library(SingleCellExperiment)
library(ggplot2)

result_path = "results/applications/hao_3_prime/"
dir.create(result_path, recursive = TRUE)

data_file = file.path(result_path, "data.RData")

if (!file.exists(data_file)) {

  hao_3_prime_data = CITEseqData::get_hao_3_prime_data("../CITEseqData/data")

  subset_proteins = function(sce, prop_nonzero) {

    pcts = Matrix::rowMeans(counts(altExp(sce)) != 0)

    altExp(sce) = altExp(sce)[pcts > prop_nonzero, ]

    return(sce)

  }

  hao_3_prime_data = subset_proteins(hao_3_prime_data, prop_nonzero = 0.5)

  hao_3_prime_data$split = paste0(hao_3_prime_data$donor, "_", hao_3_prime_data$time)
  hao_3_prime_split = lapply(unique(hao_3_prime_data$split), function(x) hao_3_prime_data[, which(hao_3_prime_data$split == x)])
  names(hao_3_prime_split) = unique(hao_3_prime_data$split)
  hao_3_prime_split = hao_3_prime_split[sort(names(hao_3_prime_split))]

  rm(hao_3_prime_data)
  gc()

  train_data = names(hao_3_prime_split)[1:12]
  validation_data = names(hao_3_prime_split)[13:18]
  test_data = names(hao_3_prime_split)[19:24]

  select_genes = function(train_sce_list, all_sce_list, n_genes) {

    genes = Reduce(intersect, lapply(all_sce_list, rownames))

    train_sce_list = lapply(train_sce_list, function(x) x[genes, ])

    ranks = sapply(train_sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

    genes = genes[order(rowMeans(ranks))][1:n_genes]

    return(genes)

  }

  genes = select_genes(hao_3_prime_split[train_data], hao_3_prime_split, 1000)

  X_list = lapply(hao_3_prime_split, function(x) t(as.matrix(logcounts(x)[genes, ])))
  Y_list = lapply(hao_3_prime_split, function(x) t(as.matrix(logcounts(altExp(x)))))
  metadata_list = lapply(hao_3_prime_split, function(x) as.data.frame(colData(x)))

  rm(hao_3_prime_split)
  gc()

  prepared_train = prepare_Y_and_indices_train(Y_list[train_data])
  prepared_validation = prepare_Y_and_indices_test(Y_list[validation_data], prepared_train$indices_list)
  prepared_test = prepare_Y_and_indices_test(Y_list[test_data], prepared_train$indices_list)

  save.image(file = data_file)

} else {

  load(data_file)

}


MultiLORS_fit = MultiLORS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                          prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list,
                          verbose = 1, n_iter = 1000, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = TRUE,
                          return_L = TRUE, n_cores = 20)

saveRDS(MultiLORS_fit, file.path(result_path, "MultiLORS_fit.rds"))

MultiLORS_refit = refit_MultiLORS(MultiLORS_fit, prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                                  prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

saveRDS(MultiLORS_refit, file.path(result_path, "MultiLORS_refit.rds"))

glmnet_fit = fit_glmnet(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                        prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

saveRDS(glmnet_fit, file.path(result_path, "glmnet_fit.rds"))

glmnet_refit = refit_glmnet(glmnet_fit, prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                            prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

saveRDS(glmnet_refit, file.path(result_path, "glmnet_refit.rds"))

OLS_fit = fit_OLS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list)

saveRDS(glmnet_refit, file.path(result_path, "OLS_fit.rds"))


(best_indices_MultiLORS = which_min(MultiLORS_fit$tuning$validation$SSE))
(best_indices_MultiLORS_refit = which_min(MultiLORS_refit$tuning$validation$SSE))
(best_index_glmnet = which_min(glmnet_fit$tuning$validation$SSE))
(best_index_glmnet_refit = which_min(glmnet_refit$tuning$validation$SSE))

MultiLORS_Beta_hat = MultiLORS_fit$model_fits[[best_indices_MultiLORS[1]]][[best_indices_MultiLORS[2]]]$Beta
MultiLORS_refit_Beta_hat = MultiLORS_refit$model_fits[[best_indices_MultiLORS_refit[1]]][[best_indices_MultiLORS_refit[2]]]$Beta
glmnet_Beta_hat = glmnet_fit$model_fits[[best_index_glmnet]]$Beta
glmnet_refit_Beta_hat = glmnet_refit$model_fits[[best_index_glmnet_refit]]$Beta
OLS_Beta_hat = OLS_fit

compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_Beta_hat)
compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, MultiLORS_refit_Beta_hat)
compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_Beta_hat)
compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, glmnet_refit_Beta_hat)
compute_avg_R2(prepared_test$Y_list, X_list[test_data], prepared_test$indices_list, prepared_train$Y_list, prepared_train$indices_list, OLS_Beta_hat)