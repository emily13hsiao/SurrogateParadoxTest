# Run A3: NNR simulations

###############################################################################
#############################      Settings      ##############################
###############################################################################
# Sample size; also 1000 and 2000
n <- 500

# Acceptable values are 1:9
setting <- 1

###############################################################################
###########################      Running Sims      ############################
###############################################################################
source("test-script-2024-09-16.R")
source("generate-data-2024-09-16.R")
set.seed(batch.num)

n.iter <- 10 # Because I am parallelizing 100 jobs

nnr_results <- list()
for (iter in 1:n.iter) {
  print(iter)
  data <- a3.data(n, setting)
  test_result <- test_assumptions(s0 = data$s0, y0 = data$y0,
                                  s1 = data$s1, y1 = data$y1,
                                  type = "nnr")
  nnr_results[[iter]] <- test_result$nnr_result
}
filename <- paste0("./a3-results/setting", toString(setting), 
                   "n", toString(n), 
                   "batch", toString(batch.num), ".RDS")
saveRDS(nnr_results, file = filename)



