# Run simulations for A2: monotonicity
library(MonotonicityTest)

# Ensure the directory exists
if (!dir.exists("./a2-results")) {
  dir.create("./a2-results")
}

###############################################################################
#############################      Settings      ##############################
###############################################################################
# Desired variance for this setting.
# Default value is 0.1, with lower variance settings they are 0.05 and 0.01
sim_sd <- 0.10

# Acceptable values are:
#   - "decreasing" = setting 6
#   - "flat" = setting 3
#   - "increasing" = setting 2
#   - "hall" = setting 4
#   - "parabola" = setting 5
#   - "increasing_linear" = setting 1

###############################################################################
###########################      Running Sims      ############################
###############################################################################
source("test-script-2025-01-20.R")
source("generate-data-2024-09-18.R")
set.seed(batch.num)

n.iter <- 10 
n <- 100

if(n==500) {ss =""}
if(n==250) {ss="250"}
if(n==100) {ss="100"}

for(set in c("increasing_linear", "parabola","hall", "increasing", "flat", "decreasing"))	{
	setting <- set	
	monotonicity_results <- list()
	for (iter in 1:n.iter) {
  		print(iter)
  		data <- a2.data(n, setting, sim_sd)
  		test_result <- test_assumptions(s0 = data$x, y0 = data$y, type = "monotonicity")
 		monotonicity_results[[iter]] <- test_result$monotonicity_result
	}
	filename <- paste0("./a2-results/", toString(setting), 
                   "sd", toString(sim_sd), 
                   "batch", ss, toString(batch.num), ".RDS")
saveRDS(monotonicity_results, file = filename)
}


