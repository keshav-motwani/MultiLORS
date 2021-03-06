library(CITEseqData)
library(SingleCellExperiment)
library(scanalysis)

result_path = "results/applications/hao_3_prime/"
dir.create(result_path, recursive = TRUE)

hao_3_prime_data = get_hao_3_prime_data("../CITEseqData/data")

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

select_genes = function(train_sce_list, all_sce_list, n_genes) {

  genes = Reduce(intersect, lapply(all_sce_list, rownames))

  train_sce_list = lapply(train_sce_list, function(x) x[genes, ])

  ranks = sapply(train_sce_list, function(x) rank(-1 * Seurat::FindVariableFeatures(assay(x, "logcounts"))$vst.variance.standardized))

  genes = genes[order(rowMeans(ranks))][1:n_genes]

  return(genes)

}

genes = select_genes(hao_3_prime_split[1:12], hao_3_prime_split, 1000)

X_list = lapply(hao_3_prime_split, function(x) t(as.matrix(logcounts(x)[genes, ])))
Y_list = lapply(hao_3_prime_split, function(x) t(as.matrix(logcounts(altExp(x)))))
metadata_list = lapply(hao_3_prime_split, function(x) as.data.frame(colData(x)))

rm(hao_3_prime_split)
gc()

saveRDS(list(X_list = X_list, Y_list = Y_list, metadata_list = metadata_list), file.path(result_path, "data.rds"))