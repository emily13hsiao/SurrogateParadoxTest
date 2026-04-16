# Sims with SE read-in and table for PAB

n.study = 10
n.each = 100

truth = read.table(paste0("truep_", n.study, n.each, ".txt"), header=T)

# Linear Table
linear_table <- matrix(NA, nrow = 6, ncol = 5)
spline <- FALSE
degree <- 1
for (setting in 1:6) {
  results <- readRDS(paste0("setting", setting, "spline", spline, "degree", degree, "results.RDS"))
  p_list <- rep(NA, 1000)
  se_list <- rep(NA, 1000)
  coverage <- rep(NA, 1000)
  for (i in 1:1000) {
    p_hat <- results[[i]]$p_hat
    se_hat <- results[[i]]$se_p
    p_list[i] <- p_hat
    se_list[i] <- se_hat
    coverage[i] <- truth[setting] >= p_hat - 1.96 * se_hat && truth[setting] <= p_hat + 1.96 * se_hat
  }
  
  linear_table[setting, ] <- c(mean(p_list), sd(p_list), mean(se_list), mean(se_list - sd(p_list)), mean(coverage))
}
linear_table

# Cubic Table
cubic_table <- matrix(NA, nrow = 6, ncol = 5)
spline <- FALSE
degree <- 3
for (setting in 1:6) {
  results <- readRDS(paste0("setting", setting, "spline", spline, "degree", degree, "results.RDS"))
  p_list <- rep(NA, 1000)
  se_list <- rep(NA, 1000)
  coverage <- rep(NA, 1000)
  for (i in 1:1000) {
    p_hat <- results[[i]]$p_hat
    se_hat <- results[[i]]$se_p
    p_list[i] <- p_hat
    se_list[i] <- se_hat
    coverage[i] <- truth[setting] >= p_hat - 1.96 * se_hat && truth[setting] <= p_hat + 1.96 * se_hat
  }
  
  cubic_table[setting, ] <- c(mean(p_list), sd(p_list), mean(se_list), mean(se_list - sd(p_list)), mean(coverage))
}
cubic_table

# Spline Table
spline_table <- matrix(NA, nrow = 6, ncol = 5)
spline <- TRUE
degree <- 3
for (setting in 1:6) {
  results <- readRDS(paste0("setting", setting, "spline", spline, "degree", degree, "results.RDS"))
  p_list <- rep(NA, 100)
  se_list <- rep(NA, 100)
  coverage <- rep(NA, 100)
  for (i in 1:100) {
    p_hat <- results[[i]]$p_hat
    se_hat <- results[[i]]$se_p
    p_list[i] <- p_hat
    se_list[i] <- se_hat
    coverage[i] <- truth[setting] >= p_hat - 1.96 * se_hat && truth[setting] <= p_hat + 1.96 * se_hat
  }
  
  spline_table[setting, ] <- c(mean(p_list), sd(p_list), mean(se_list), mean(se_list - sd(p_list)), mean(coverage))
}
spline_table

# Full table
full_table <- round(cbind(truth, linear_table, cubic_table, spline_table), 3)
full_table





