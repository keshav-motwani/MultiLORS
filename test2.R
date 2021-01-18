library(MultiLORS)
library(CITEseqData)
library(SingleCellExperiment)

datasets = get_all_datasets(cache_path = "../CITEseqData/data/")

train_data = c("kotliarov_data", "stoeckius_cbmc_data", "pbmc_1k_protein_v3_data", "malt_10k_protein_v3_data")
validation_data = c("5k_pbmc_protein_v3_nextgem_data", "pbmc_10k_protein_v3_data")
test_data = c("SC3_v3_NextGem_DI_PBMC_CSP_1K_data", "stoeckius_pbmc_data")

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

prepared_train = prepare_Y_and_indices_train(Y_list[train_data])
prepared_validation = prepare_Y_and_indices_test(Y_list[validation_data], prepared_train$indices_list)
prepared_test = prepare_Y_and_indices_test(Y_list[test_data], prepared_train$indices_list)

fit_1 = fit_glmnet(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                   prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

fit_2 = MultiLORS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                  prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list,
                  verbose = 1, n_iter = 100, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = TRUE)
