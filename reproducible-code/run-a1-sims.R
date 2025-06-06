# Run simulations for A1: stochastic dominance


# Ensure the directory exists
if (!dir.exists("./a1-results")) {
  dir.create("./a1-results")
}

###############################################################################
#############################      Settings      ##############################
###############################################################################
n <- 500
n.iter <- 1000

###############################################################################
###########################      Running Sims      ############################
###############################################################################
source("test-script-2024-09-18.R")
source("generate-data-2024-09-18.R")
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



