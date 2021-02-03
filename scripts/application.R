library(MultiLORS)
library(CITEseqData)
library(SingleCellExperiment)

datasets = get_all_datasets(cache_path = "../CITEseqData/data/")

subset_proteins = function(sce, prop_nonzero) {

  pcts = Matrix::rowMeans(counts(altExp(sce)) != 0)

  altExp(sce) = altExp(sce)[pcts > prop_nonzero, ]

  return(sce)

}

datasets = lapply(datasets, subset_proteins, prop_nonzero = 0.5)

datasets$hao_3_prime_data$split = paste0(datasets$hao_3_prime_data$donor, "_", datasets$hao_3_prime_data$time)
hao_3_prime_split = lapply(unique(datasets$hao_3_prime_data$split), function(x) datasets$hao_3_prime_data[, which(datasets$hao_3_prime_data$split == x)])
names(hao_3_prime_split) = paste0("hao_3_prime_data_", unique(datasets$hao_3_prime_data$split))
datasets$hao_3_prime_data = NULL

datasets$hao_5_prime_data$split = paste0(datasets$hao_5_prime_data$donor, "_", datasets$hao_5_prime_data$time)
hao_5_prime_split = lapply(unique(datasets$hao_5_prime_data$split), function(x) datasets$hao_5_prime_data[, which(datasets$hao_5_prime_data$split == x)])
names(hao_5_prime_split) = paste0("hao_5_prime_data_", unique(datasets$hao_5_prime_data$split))
datasets$hao_5_prime_data = NULL

datasets = c(datasets, hao_3_prime_split, hao_5_prime_split)

validation_data = "su_data"
test_data = "kotliarov_data"
train_data = setdiff(names(datasets), c(validation_data, test_data))

select_genes = function(train_sce_list, all_sce_list, n_genes) {

  genes = Reduce(intersect, lapply(all_sce_list, rownames))

  train_sce_list = lapply(train_sce_list, function(x) x[genes, ])

  ranks = sapply(train_sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))][1:n_genes]

  return(genes)

}

genes = select_genes(datasets[train_data], datasets, 1000)

X_list = lapply(datasets, function(x) t(as.matrix(logcounts(x)[genes, ])))
Y_list = lapply(datasets, function(x) t(as.matrix(logcounts(altExp(x)))))

rm(datasets)
gc()

prepared_train = prepare_Y_and_indices_train(Y_list[train_data])
prepared_validation = prepare_Y_and_indices_test(Y_list[validation_data], prepared_train$indices_list)
prepared_test = prepare_Y_and_indices_test(Y_list[test_data], prepared_train$indices_list)

fit_2 = MultiLORS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                  prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list,
                  verbose = 1, n_iter = 500, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = TRUE,
                  return_L = FALSE, n_cores = 20)

saveRDS(fit_2, "MultiLORS_fit.rds")

fit_1 = fit_glmnet(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                   prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

saveRDS(fit_1, "glmnet_fit.rds")

