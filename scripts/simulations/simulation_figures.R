library(tidyverse)
library(patchwork)

RESULT_PATH = "results/simulations_new"

files = list.files(RESULT_PATH, full.names = TRUE)
files = files[grepl("rds", files)]

results = lapply(files, function(x) {
  result = readRDS(x)
  data.frame(experiment = result$parameters$experiment,
             value = result$parameters[[result$parameters$experiment]],
             method = result$parameters$method,
             replicate = result$parameters$replicate,
             Beta_SSE = result$result$Beta_SSE,
             test_R2 = result$result$test_R2,
             test_SSE = result$result$test_SSE,
             test_SST = result$result$test_SST)
})

result = do.call(rbind, results)

methods = c("MultiLORS", "glmnet", "ORC_L_glmnet", "ORC_ALL_MultiLORS", "ORC_ALL_glmnet", "ORC_L_ALL_glmnet")
methods = c(methods, setdiff(unique(result$method), methods))

summary = result %>%
  pivot_longer(Beta_SSE:test_SST, names_repair = "minimal", values_to = "result") %>%
  mutate(value = factor(value)) %>%
  mutate(method = factor(method, levels = methods)) %>%
  group_by(experiment, value, method, name) %>%
  summarize(mean = mean(result), two_se = 2 * sd(result)/sqrt(n())) %>%
  mutate(group = paste0(experiment, method, name))

plasma_pal = viridis::plasma(n = length(methods) + 2)[1:length(methods)]

plots = list()

for (exp in unique(summary$experiment)) {
  plots = c(
    plots,
    list(
      ggplot(
        summary %>% filter(experiment == exp),
        aes(
          x = value,
          y = mean,
          color = method,
          group = method,
          linetype = method,
          ymin = mean - two_se,
          ymax = mean + two_se
        )
      ) +
        geom_point(size = 0.5, position=position_dodge(0.2)) +
        geom_line(position=position_dodge(0.2)) +
        geom_errorbar(width = 0.2, position=position_dodge(0.2), linetype = "solid", show.legend = FALSE) +
        facet_wrap( ~ name, scales = "free", nrow = 1) +
        theme_classic() +
        theme(strip.background = element_blank(), strip.placement = "outside") +
        scale_color_manual(values = plasma_pal) +
        theme(legend.position = "bottom") +
        xlab(exp) +
        ylab(NULL)
    )
  )
}
names(plots) = unique(summary$experiment)

plots[["r"]] =plots[["r"]] + xlab("Rank of L")
plots[["R2"]] = plots[["R2"]] + xlab(expression(Population~R^{2}))
plots[["theta"]] = plots[["theta"]] + xlab( expression(theta))

wrap_plots(plots, nrow = length(plots)) + plot_layout(guides = "collect") &
  theme(legend.position='bottom')

dir.create(file.path(RESULT_PATH, "figures"))
ggsave(file.path(RESULT_PATH, "figures", "simulation_figures.pdf"), height = 4.2 * length(unique(result$experiment)), width = 16)
write.csv(result, file.path(RESULT_PATH, "figures", "simulation_results.csv"))
write.csv(summary, file.path(RESULT_PATH, "figures", "simulation_result_summary.csv"))
