# Run simulations for A1: stochastic dominance

###############################################################################
#############################      Settings      ##############################
###############################################################################
n <- 2000
n.iter <- 1000

###############################################################################
###########################      Running Sims      ############################
###############################################################################
source("test-script-2024-09-16.R")
source("generate-data-2024-09-16.R")
set.seed(1)

for (setting in 1:9) {
  sd_results <- list()
  for (iter in 1:n.iter) {
    print(iter)
    data <- a1.data(n, setting)
    test_result <- test_assumptions(s0 = data$s0, s1 = data$s1, type = "sd")
    sd_results[[iter]] <- test_result$sd_result
  }
  filename <- paste0("./a1-results/setting",
                     toString(setting), "n", toString(n), ".RDS")
  saveRDS(sd_results, file = filename)
}



