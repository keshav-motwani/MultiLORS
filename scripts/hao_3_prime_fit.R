library(MultiLORS)
library(CITEseqData)
library(SingleCellExperiment)
library(ggplot2)

result_path = "results/applications/hao_3_prime/"
dir.create(result_path, recursive = TRUE)

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

rm(hao_3_prime_split)
gc()

prepared_train = prepare_Y_and_indices_train(Y_list[train_data])
prepared_validation = prepare_Y_and_indices_test(Y_list[validation_data], prepared_train$indices_list)
prepared_test = prepare_Y_and_indices_test(Y_list[test_data], prepared_train$indices_list)

save.image(file = file.path(result_path, "data.RData"))

MultiLORS_fit = MultiLORS(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                  prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list,
                  verbose = 1, n_iter = 500, tolerance = 1e-6, n_lambda = 20, n_gamma = 20, early_stopping = TRUE,
                  return_L = FALSE, n_cores = 20)

saveRDS(MultiLORS_fit, file.path(result_path, "MultiLORS_fit.rds"))

glmnet_fit = fit_glmnet(prepared_train$Y_list, X_list[train_data], prepared_train$indices_list,
                   prepared_validation$Y_list, X_list[validation_data], prepared_validation$indices_list)

saveRDS(glmnet_fit, file.path(result_path, "glmnet_fit.rds"))