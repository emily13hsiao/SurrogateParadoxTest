# Objective: read in simulation files

full_results <- data.frame()
for (batch.num in 1:20) {
  filename <- paste0("./run-simulations/polynomial/setting6/batch", batch.num, ".RDS")
  results <- data.frame(readRDS(filename))
  full_results <- rbind(full_results, results)
}

df <- data.frame(full_results)
colnames(df) <- c("p_hat", "p_hat_se", "q_hat", "q_hat_se", "true_DeltaB")

# True Delta B
mean(df$true_DeltaB < 0) # 0.631
mean(df$p_hat) # 0.56

# Empirical SD of p_hat compared to estimated values
sd(df$p_hat) # 0.29
mean(df$p_hat_se) # 0.245

# Mean q_hat
quantile(df$true_DeltaB, 0.1)[[1]] # -1.2
mean(df$q_hat) # -1.03

# Empirical SD of q_hat
sd(df$q_hat) # 0.7255
mean(df$q_hat_se) # 0.712






