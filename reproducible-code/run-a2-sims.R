# Run simulations for A2: monotonicity

###############################################################################
#############################      Settings      ##############################
###############################################################################
# Desired variance for this setting.
# Default value is 0.1, with lower variance settings they are 0.05 and 0.01
sim_sd <- 0.01

# Acceptable values are:
#   - "decreasing"
#   - "flat"
#   - "increasing"
#   - "hall"
#   - "parabola"
#   - "increasing_linear"
setting <- "increasing_linear"

###############################################################################
###########################      Running Sims      ############################
###############################################################################
source("test-script-2024-09-16.R")
source("generate-data-2024-09-18.R")
set.seed(batch.num)

n.iter <- 10 # Because I am parallelizing 100 jobs
n <- 500

monotonicity_results <- list()
for (iter in 1:n.iter) {
  print(iter)
  data <- a2.data(n, setting, sim_sd)
  test_result <- test_assumptions(s0 = data$x, y0 = data$y, type = "monotonicity")
  monotonicity_results[[iter]] <- test_result$monotonicity_result
}
filename <- paste0("./a2-results/", toString(setting), 
                   "sd", toString(sim_sd), 
                   "batch", toString(batch.num), ".RDS")
saveRDS(monotonicity_results, file = filename)



