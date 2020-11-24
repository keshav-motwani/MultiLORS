library(MultiLORS)

ARRAY_ID = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
RESULT_PATH = "results/simulations"
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
  R2 = seq(from = 0.5, to = 0.9, by = 0.05),
  r = 2:15,
  theta = seq(from = 10, to = 100, by = 10)
)

parameters = expand_parameters(considered_values, defaults, 50, c("MultiLORS", "glmnet", "ORC_L_glmnet", "ORC_ALL_MultiLORS", "ORC_L_ALL_glmnet"))[[ARRAY_ID]]
print(parameters)
result = evaluate_parameters(parameters)
saveRDS(result, file.path(RESULT_PATH, paste0(parameters$experiment, "_", gsub(".", "_", parameters[[parameters$experiment]], fixed = TRUE), "_", parameters$method, "_", parameters$replicate, ".rds")))
print(result$result$test_R2)
print(result$result$Beta_SSE)