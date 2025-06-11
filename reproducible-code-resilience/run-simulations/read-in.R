# Objective: read in simulation files

setting <- 8

full_results <- data.frame()
for (batch.num in 1:20) {
  filename <- paste0("setting",setting, "batch", batch.num, ".RDS")
  results <- data.frame(readRDS(filename))
  full_results <- rbind(full_results, results)
}

df <- data.frame(full_results)
colnames(df) <- c("p_hat", "p_hat_se", "q_hat", "q_hat_se", "true_DeltaB")

# True Delta B
mean(df$true_DeltaB < 0)
mean(df$p_hat) 

# Empirical SD of p_hat compared to estimated values
sd(df$p_hat) 
mean(df$p_hat_se) 

# Mean q_hat
quantile(df$true_DeltaB, 0.1)[[1]]
mean(df$q_hat) 

# Empirical SD of q_hat
sd(df$q_hat) 
mean(df$q_hat_se) 






