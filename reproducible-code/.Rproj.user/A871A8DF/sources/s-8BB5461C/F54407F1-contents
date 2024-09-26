# Read in a1 sims

n <- 500

###############################################################################
###########################      Read Results      ############################
###############################################################################

for (setting in 1:9) {
  results <- readRDS(paste0("./a1-results/setting", toString(setting), 
                            "n", toString(n), ".RDS"))
  p_vals <- sapply(1:n.iter, function(ii) results[[ii]]$p.value)
  p_sd <- mean(p_vals > 0.05)
  print(paste("Setting", toString(setting), "P_SD:", toString(p_sd)))
}


