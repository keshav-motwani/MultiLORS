library(MultiLORS)

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/simulations_new"
dir.create(RESULT_PATH, recursive = TRUE)

defaults = list(
  R2 = 0.75,
  r = 5,
  theta = 50,
  sigma_2 = 1,
  sparsity = 0.9,
  p = 500,
  q = 50,
  O_k_train = ceiling(seq(from = 0.7, to = 0.2, length.out = 5) * 50),
  n_k_train = seq(from = 100, to = 500, length.out = 5),
  O_k_validation = c(25, 25),
  n_k_validation = c(300, 300),
  O_k_test = rep(50, 5),
  n_k_test = rep(1000, 5)
)

considered_values = list(
  r = 2:15,
  R2 = seq(from = 0.5, to = 0.9, by = 0.05),
  theta = seq(from = 10, to = 100, by = 10)
)

methods = c("MultiLORS", "glmnet", "ORC_L_MultiLORS", "ORC_L_glmnet", "ORC_ALL_MultiLORS", "ORC_ALL_glmnet", "ORC_L_ALL_MultiLORS", "ORC_L_ALL_glmnet")

parameters = expand_parameters(considered_values, defaults, 50, methods)

for (i in 1:length(methods)) {

  current_parameters = parameters[[ARRAY_ID * length(methods) + i]]

  result = evaluate_parameters(current_parameters)
  saveRDS(result, file.path(RESULT_PATH, paste0(current_parameters$experiment, "_", gsub(".", "_", current_parameters[[current_parameters$experiment]], fixed = TRUE), "_", current_parameters$method, "_", current_parameters$replicate, ".rds")))

}