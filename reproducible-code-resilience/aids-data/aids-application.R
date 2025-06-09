# Objective: Run methods on the aids data.
library(tidyverse)
library(MASS)
library(patchwork)
source("./../functions.R")
set.seed(1)

# Read in data
s0.A <- readRDS("./s0_A.RDS")
y0.A <- readRDS("./y0_A.RDS")
s1.A <- readRDS("./s1_A.RDS")
y1.A <- readRDS("./y1_A.RDS")
s0.B <- readRDS("./s0_B.RDS")
s1.B <- readRDS("./s1_B.RDS")

Delta_A <- mean(y1.A) - mean(y0.A)
Delta_A

################################################################################
############################ Resilience Intervals ##############################
################################################################################

# GP FIXED
gp_fixed_results <- gaussian_process_interval(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, 
                                              sigma2 = 1.25, theta = 2, plot = TRUE)
gp_fixed_results$p_hat
gp_fixed_results$p_se
gp_fixed_results$q_hat
gp_fixed_results$q_se

# FOURIER
# So first let me calculate the average residuals which is what I am using to choose
# the variance parameters.
mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))
average_residual <- mean(c(abs(sapply(s0.A, mu_hat_0) - y0.A), abs(sapply(s1.A, mu_hat_1) - y1.A)))
variance_vector <- average_residual * c(0.5, 0.5, 0.25, 0.25)
period_vector <- c(0.5, 0.25, 0.1)
print(variance_vector)
fourier_results <- fourier_interval(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B,
                                    var_vec = variance_vector,
                                    period = period_vector, plot = TRUE)
fourier_results$p_hat
fourier_results$p_se
fourier_results$q_hat
fourier_results$q_se

# POLYNOMIAL
variance_vector <- average_residual * c(0.1, 0.1, 0.025, 0.025)
print(variance_vector)
polynomial_results <- polynomial_interval(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B,
                                          var_vec = variance_vector, plot = TRUE)
polynomial_results$p_hat
polynomial_results$p_se
polynomial_results$q_hat
polynomial_results$q_se

################################################################################
############################### Resilience Sets ################################
################################################################################

sigma2_vals <- seq(0.1, 10, length.out = 50); theta_vals <- seq(0.1, 10, length.out = 50)
gp_set <- gp_resilience_set(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sigma2_vals, theta_vals)
gp_set
ggsave("aids_gp_res_set.png")

sig1_vals <- seq(0.01, 0.25, length.out = 50); sig2_vals <- seq(0.01, 2, length.out = 50)
fourier_set <- fourier_resilience_set(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sig1_vals, sig2_vals)
fourier_set
ggsave("aids_fourier_res_set.png")

sig1_vals <- seq(0.01, 0.15, length.out = 50); sig2_vals <- seq(0.01, 0.2, length.out = 50)
poly_set <- polynomial_resilience_set(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sig1_vals, sig2_vals)
poly_set
ggsave("aids_poly_res_set.png")

################################################################################
################################# Data Plots ###################################
################################################################################

# GP VERSION
control_plot <- gp_fixed_results$control_plot
control_plot <- control_plot + labs(title = "Control Group",
                                    x = "Change in CD4 count",
                                    y = "Change in RNA HIV-1 concentration")
treatment_plot <- gp_fixed_results$treatment_plot
treatment_plot <- treatment_plot + labs(title = "Treatment Group",
                                        x = "Change in CD4 count",
                                        y = "Change in RNA HIV-1 concentration")
treatment_plot + control_plot + plot_layout(guides = "collect") &
  theme(axis.title.x = element_text(), axis.title.y = element_text())
ggsave("aids_gp_plot.png")

# FOURIER VERSION

control_plot <- fourier_results$control_plot
control_plot <- control_plot + labs(title = "Control Group",
                                    x = "Change in CD4 count",
                                    y = "Change in RNA HIV-1 concentration")
treatment_plot <- fourier_results$treatment_plot
treatment_plot <- treatment_plot + labs(title = "Treatment Group",
                                        x = "Change in CD4 count",
                                        y = "Change in RNA HIV-1 concentration")
treatment_plot + control_plot + plot_layout(guides = "collect") &
  theme(axis.title.x = element_text(), axis.title.y = element_text())
ggsave("aids_fourier_plot.png")

# POLYNOMIAL VERSION

control_plot <- polynomial_results$control_plot
control_plot <- control_plot + labs(title = "Control Group",
                                    x = "Change in CD4 count",
                                    y = "Change in RNA HIV-1 concentration")
treatment_plot <- polynomial_results$treatment_plot
treatment_plot <- treatment_plot + labs(title = "Treatment Group",
                                        x = "Change in CD4 count",
                                        y = "Change in RNA HIV-1 concentration")
treatment_plot + control_plot + plot_layout(guides = "collect") &
  theme(axis.title.x = element_text(), axis.title.y = element_text())
ggsave("aids_polynomial_plot.png")









