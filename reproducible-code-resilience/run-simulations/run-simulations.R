# Objective: run simulations
library(tidyverse)
library(SurrogateParadoxTest)
library(MASS)
source("./../functions.R")

setting <- 1 # Choices are: 1-10
method <- "gp_fixed" # gp_fixed, polynomial, fourier
sim_reps <- 50
n.A <- 400
n.B <- 200
get_var <- TRUE

# Additional Settings
if (method == "gp_fixed" | method == "gp_estimated") {
  if (setting == 1) {
    sigma2 <- 0.3
    theta <- 5
  } else if (setting == 2) {
    sigma2 <- 0.25
    theta <- 5
  } else if (setting == 3) {
    sigma2 <- 1
    theta <- 1
  }
} else if (method == "polynomial") {
  if (setting == 4) {
    var_vec = c(0.25, 0.25, 0.1, 0.1)
  } else if (setting == 5) {
    var_vec = c(0.25, 0.25, 0.1, 0.1)
  } else if (setting == 6) {
    var_vec = c(0.25, 0.25, 0.1, 0.1)
  }
} else if (method == "fourier") {
  if (setting == 7) {
    period = c(0.5, 0.25, 0.1)
    var_vec = c(0.5, 0.5, 0.1, 0.1)
  } else if (setting == 8) {
    period = c(0.5, 0.25, 0.1)
    var_vec = c(0.5, 0.5, 0.1, 0.1)/10
  } else if (setting == 9) {
    period = c(0.5, 0.25, 0.1)
    var_vec = c(0.5, 0.5, 0.1, 0.1)/10
  }
}

set.seed(batch.num)

# Table to store results
all_data_table <- matrix(data = NA, nrow = sim_reps, ncol = 5)
colnames(all_data_table) <- c("p_hat", "p_se", "q_hat", "q_se", "true_DeltaB")
for (jj in 1:sim_reps) {

  print(c("iter", jj))

  # Generate Study B
  if (method == "gp_fixed") {
    data <- generate_data(setting, n.A, n.B)
  } else if (method == "polynomial") {
    data <- generate_data(setting, n.A, n.B, var_vec)
  } else if (method == "fourier") {
    data <- generate_data(setting, n.A, n.B, var_vec, period)
  }

  # Run the results
  if (method == "gp_fixed") {
    result <- gaussian_process_interval(data$s0.A, data$y0.A, data$s1.A, data$y1.A,
                                      data$s0.B, data$s1.B, sigma2 = sigma2, theta = theta)
  }else if (method == "polynomial") {
    result <- polynomial_interval(data$s0.A, data$y0.A, data$s1.A, data$y1.A,
                                  data$s0.B, data$s1.B, var_vec)
  } else if (method == "fourier") {
    result <- fourier_interval(data$s0.A, data$y0.A, data$s1.A, data$y1.A, 
                                       data$s0.B, data$s1.B, var_vec, period)
  }
  
  # Extract the numbers from the results
  p_hat <-result$p_hat
  p_se <- result$p_se
  q_hat <- result$q_hat
  q_se <- result$q_se
  true_DeltaB <- data$Delta_B

  # Add these to the table
  row_vec <- c(p_hat, p_se, q_hat, q_se, true_DeltaB)
  all_data_table[jj, ] <- row_vec
}

saveRDS(all_data_table, file = paste0("./all_results/setting", setting, "/batch", batch.num, ".RDS"))




