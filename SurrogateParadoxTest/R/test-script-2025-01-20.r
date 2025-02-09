# File for all the tests

###############################################################################
###########################      Assumption 1      ############################
###############################################################################

barrett_donald_cutoff <- function(alpha) {
  sqrt(-1 * log(alpha) / 2)
}

barrett_donald_p <- function(statistic) {
  exp(-2 * statistic^2)
}

sd_test <- function(s0, s1, alpha = 0.05) {
  # Make grid
  combined_sample <- c(s0, s1)
  min.x <- min(combined_sample)
  max.x <- max(combined_sample)
  grid <- seq(min.x, max.x, length.out = length(s0) + length(s1))
  
  # Calculate ECDFs
  ecdf0 <- ecdf(s0)
  ecdf1 <- ecdf(s1)
  ecdf_diff <- function(x) ecdf1(x) - ecdf0(x)
  
  # Calculate the max difference
  grid.y <- sapply(grid, ecdf_diff)
  sup <- max(grid.y)
  s_hat <- sqrt( (length(s0) * length(s1)) / (length(s0) + length(s1))) * sup
  
  # Test decision
  cutoff <- barrett_donald_cutoff(alpha)
  p.value <- barrett_donald_p(s_hat)
  reject <- p.value <= alpha
  
  return(list(s_hat = s_hat, p.value = p.value, reject = reject))
}

###############################################################################
#########################      Helper Functions      ##########################
###############################################################################

gaussian_kernel <- function(x) {
  exp(-x^2 / 2) / sqrt(2 * pi)
}

# Takes in a dataset and a kernel function
# Outputs a function which is the kernel smoother
# Parameter to this function is a vector of regressors
smoother_fitter <- function(X, Y, kernel = gaussian_kernel, h) {
  smoother <- function(x) {
    weights = sapply(X, function(val) kernel((x - val)/h) /h )
    sum(weights * Y) / sum(weights)
  }
  smoother
}

calculate_bandwidth <- function(s) {
  bw.nrd(s) * length(s)^(-0.1)
}

###############################################################################
###########################      Assumption 2      ############################
###############################################################################

S <- function(a, b, r, s, X, Y) {
  s <- 0
  for (i in r:s) {
    s <- s + (Y[i] - (a + b * X[i]))^2
  }
  s
}

Q <- function(r, s, X) {
  subset_X <- X[r:s]
  avg_X <- sum(subset_X) / (s - r + 1)
  sqrt(sum((subset_X - avg_X)^2))
}

a_b_hat <- function(r, s, X, Y) {
  subset_X <- X[r:s]
  subset_Y <- Y[r:s]
  df <- data.frame(X = subset_X, Y = subset_Y)
  coeffs <- lm(Y ~ X, data = df)$coefficients
  a_hat <- coeffs[1]; b_hat <- coeffs[2]
  return(list(a_hat = a_hat, b_hat = b_hat))
}

T_m <- function(m, X, Y) {
  n <- length(X)
  stat_vals <- c()
  b_vals <-c()
  Q_vals <-c()
  for (r in 1:(n - m + 1)) { 
    for (s in (r + m - 1):n) {
      b_hat <- a_b_hat(r, s, X, Y)$b_hat
      Q_r_s <- Q(r, s, X)
      b_vals = c(b_vals, b_hat)
      Q_vals <- c(Q_vals, Q_r_s)
      stat_vals <- append(stat_vals, -1 * b_hat * Q_r_s)
    }
  }
  return(list("stat" = max(stat_vals, na.rm = T), "stat_vals" = stat_vals, 
              "b_vals" = b_vals, "Q_vals" = Q_vals)
  )
}

monotonicity_test <- function(X, Y, h = NA, m = 5, bootstrap_n = 100, 
                              alpha = 0.05) {
  # Calculate actual T_m
  T_m_actual <- T_m(m, X, Y)
  
  # Bootstrap it
  n <- length(X)
  g_hat <- smoother_fitter(X, Y, h = bw.nrd(X)*(n^(-0.1)))
  preds <- sapply(1:n, function(i) g_hat(X[i]))
  errors <- sapply(1:n, function(i) Y[i] - g_hat(X[i]))
  T_m_samples <- rep(NA, bootstrap_n)
  
    for (i in 1:bootstrap_n) {
      resampled_ind <- sample(1:n, size = n, replace = TRUE)
      resampled_errors <- errors[resampled_ind]
      resampled_X = X[resampled_ind]
      
      new.x = resampled_X[order(resampled_X)]
      new.y = resampled_errors[order(resampled_X)]
      T_m_star <- T_m(m, new.x, new.y)
      T_m_samples[i] <- T_m_star$stat
    }
  
  
  cutoff <- quantile(T_m_samples, 1 - alpha)
  
  p_val <- sum(T_m_samples > T_m_actual$stat) / bootstrap_n
  reject <- T_m_actual$stat > cutoff
  
  return(list(T_m_value = T_m_actual$stat, p_val = p_val, reject = reject, 
              T_m_samples = T_m_samples))
}

###############################################################################
###########################      Assumption 3      ############################
###############################################################################

modified_S_stat <- function(mu0_hat, mu1_hat, s0, y0, s1, y1,
                            grid_x, boot = FALSE) {
  
  h0 <- calculate_bandwidth(s0)
  h1 <- calculate_bandwidth(s1)
  
  if (boot == TRUE) {
    # Smoothed functions
    smoothed_control <- smoother_fitter(s0, y0, gaussian_kernel, h0)
    smoothed_treatment <- smoother_fitter(s1, y1, gaussian_kernel, h1)
  } else { # If not boot strapping, these should be 0
    smoothed_control <- function(x) 0
    smoothed_treatment <- function(x) 0
  }
  
  # Find max difference on the support
  difference_function <- function(x) {
    mu0_hat(x) - mu1_hat(x) - (smoothed_control(x) - smoothed_treatment(x))
  }
  
  # Find max value
  grid_y <- sapply(grid_x, difference_function)
  sup <- max(grid_y)
  
  # Calculate statistic
  s_hat <- sqrt( (length(s0) * length(s1)) / (length(s0) + length(s1))) * sup
  
  return(list(s_hat = s_hat, sup = sup))
}

nnr_test <- function(s0, y0, s1, y1, n_bootstrap = 200, alpha = 0.05) {
  # Create grid on support
  combined_x <- c(s0, s1)
  min_x <- min(combined_x)
  max_x <- max(combined_x)
  grid_x <- seq(min_x, max_x, length.out = length(combined_x))
  
  # Calculate mu0 and mu1 estimates
  h0 <- calculate_bandwidth(s0)
  h1 <- calculate_bandwidth(s1)
  mu0_hat <- smoother_fitter(s0, y0, gaussian_kernel, h0)
  mu1_hat <- smoother_fitter(s1, y1, gaussian_kernel, h1)
  s_hat <- modified_S_stat(mu0_hat, mu1_hat, s0, y0, s1, y1, grid_x)$s_hat
  
  S_vec <- rep(NA, n_bootstrap)
  
  # Bootstrap
  for (i in 1:n_bootstrap) {
    # Should this be true or false? 
    control_index <- sample(1:length(s0), length(s0), replace = TRUE)
    S_control_sample <- s0[control_index]
    Y_control_sample <- y0[control_index]
    
    treatment_index <- sample(1:length(s1), length(s1), replace = TRUE)
    S_treatment_sample <- s1[treatment_index]
    Y_treatment_sample <- y1[treatment_index]
    
    # Calculate new S stat
    new_S_hat <- modified_S_stat(mu0_hat, mu1_hat, 
                                 S_control_sample, Y_control_sample,
                                 S_treatment_sample, Y_treatment_sample,
                                 grid_x, boot = TRUE)$s_hat
    S_vec[i] <- new_S_hat
  }
  
  # Calculate p-value
  prop <- sum(S_vec >= s_hat) / length(S_vec)
  reject <- prop <= alpha
  
  # Return results
  return(list(p_value = prop, reject = reject, s_hat = s_hat, s_vec = S_vec))
}

# maybe a parameter to control the bootstrap sample size
test_assumptions <- function(s0 = NULL, y0 = NULL, s1 = NULL, y1 = NULL, 
                             trim = 0.95, alpha = 0.05, type = "all", 
                             all_results = TRUE, direction = "positive",
                             monotonicity_bootstrap_n = 100,
                             nnr_bootstrap_n = 200) {
  
  # First put the data in order - this is important for monotonicity test
  if (!is.null(s0)) {
    order0 <- order(s0)
    s0 <- s0[order0]
    if (!is.null(y0)) {
      y0 <- y0[order0]
    }
  }
  if (!is.null(s1)) {
    order1 <- order(s1)
    s1 <- s1[order1]
    if (!is.null(y1)) {
      y1 <- y1[order1]
    }
  }
  
  # Then trim the data
  if (!is.null(s0)) {
    keep_s0 <- s0 > quantile(s0, (1 - trim) / 2) & s0 < quantile(s0, 1 - (1 - trim) / 2)
    s0 <- s0[keep_s0]
    if (!is.null(y0)) {
      y0 <- y0[keep_s0]
    }
  }
  if (!is.null(s1)) {
    keep_s1 <- s1 > quantile(s1, (1 - trim) / 2) & s1 < quantile(s1, 1 - (1 - trim) / 2)
    s1 <- s1[keep_s1]
    if (!is.null(y1)) {
      y1 <- y1[keep_s1]
    }
  }
  
  # A1: stochastic dominance
  if (type == "sd" | type == "all") {
    if (direction == "positive") {
      sd_result <- sd_test(s0, s1, alpha = alpha)
    } else {
      sd_result <- sd_test(s1, s0, alpha = alpha)
    }
    
  }
  
  # A2: monotonicity
  if (type == "monotonicity" | type == "all") {
    if (!is.null(s0) & !is.null(y0)) {
    monotonicity0_result = 	c()
      if (direction == "positive") {
        temp <- MonotonicityTest::monotonicity_test(s0, y0, m = floor(0.05 * length(s0)), boot_num = monotonicity_bootstrap_n)
      } else {
        temp <- MonotonicityTest::monotonicity_test(s0, -1 * y0, m = floor(0.05 * length(s0)), boot_num = monotonicity_bootstrap_n)
      }
      monotonicity0_result$T_m_value = temp$stat
      monotonicity0_result$p_val = temp$p
      monotonicity0_result$reject = (temp$p<alpha)
      monotonicity0_result$T_m_samples = temp$dist
              
    }
    if (!is.null(s1) & !is.null(y1)) {
    	monotonicity1_result = 	c()
      if (direction == "positive") {
        temp <- MonotonicityTest::monotonicity_test(s1, y1, m = floor(0.05 * length(s1)), boot_num = monotonicity_bootstrap_n)
      } else {
        temp <- MonotonicityTest::monotonicity_test(s1, -1 * y1, m = floor(0.05 * length(s1)), boot_num = monotonicity_bootstrap_n)
      }
      monotonicity1_result$T_m_value = temp$stat
      monotonicity1_result$p_val = temp$p
      monotonicity1_result$reject = (temp$p<alpha)
      monotonicity1_result$T_m_samples = temp$dist
 
    }
  }
  
  # A3: non-negative residual treatment effect
  if (type == "nnr" | type == "all") {
    if (direction == "positive") {
      nnr_result <- nnr_test(s0, y0, s1, y1, alpha = alpha,
                             n_bootstrap = nnr_bootstrap_n)
    } else {
      nnr_result <- nnr_test(s1, y1, s0, y0, alpha = alpha,
                             n_bootstrap = nnr_bootstrap_n)
    }
  }
  
  # Make table 
  if (type == "all") {
    names_vec <- c("Stochastic dominance assumption", 
                   "Monotonicity assumption (control)",
                   "Monotonicity assumption (treatment)", 
                   "Non-negative residual treatment effect")
    results_vec <- c(sd_result$reject, monotonicity0_result$reject, 
                     monotonicity1_result$reject, nnr_result$reject)
    results_vec[results_vec == TRUE] <- "Does not hold"
    results_vec[results_vec == FALSE] <- "Holds"
    result <- cbind(names_vec, results_vec)
    colnames(result) <- c("Assumption", "Result")
    rownames(result) <- 1:4
  } else {
    if (type == "sd") {
      if (sd_result$reject) {
        result <- "Stochastic dominance assumption: Does not hold"
      } else {
        result <- "Stochastic dominance assumption: Holds"
      }
    } else if (type == "monotonicity") {
      if (!is.null(s0) & !is.null(s1)) {
        names_vec <- c("Monotonicity assumption (control)", 
                       "Monotonicity assumption (treatment)"
        )
        results_vec <- c(monotonicity0_result$reject, monotonicity1_result$reject)
        results_vec[results_vec == TRUE] <- "Does not hold"
        results_vec[results_vec == FALSE] <- "Holds"
        result <- cbind(names_vec, results_vec)
        colnames(result) <- c("Assumption", "Result")
        rownames(result) <- 1:2
      } else if (!is.null(s0)) {
        if (monotonicity0_result$reject) {
          result <- "Monotonicity assumption: Does not hold"
        } else {
          result <- "Monotonicity assumption: Holds"
        }
      } else if (!is.null(s1)) {
        if (monotonicity1_result$reject) {
          result <- "Monotonicity assumption: Does not hold"
        } else {
          result <- "Monotonicity assumption: Holds"
        }
      }
      
    } else if (type == "nnr") {
      if (nnr_result$reject) {
        result <- "Non-negative residual treatment effect assumption: Does not hold"
      } else {
        result <- "Non-negative residual treatment effect assumption: Holds"
      }
    }
  }
  
  # Return all results
  if (all_results) {
    if (type == "all") {
      return(list(
        result = result,
        sd_result = sd_result,
        monotonicity0_result = monotonicity0_result,
        monotonicity1_result = monotonicity1_result,
        nnr_result = nnr_result
      ))
    } else if (type == "sd") {
      return(list(result = result, sd_result = sd_result))
    } else if (type == "monotonicity") {
      if (!is.null(s0) & !is.null(s1)) {
        return(list(result = result, monotonicity0 = monotonicity0_result,
                    monotonicity1 = monotonicity1_result))
      } else if (!is.null(s0)) {
        return(list(result = result, monotonicity_result = monotonicity0_result))
      } else if (!is.null(s1)) {
        return(list(result = result, monotonicity_result = monotonicity1_result))
      }
    } else if (type == "nnr") {
      return(list(result = result, nnr_result = nnr_result))
    }
  } else {
    return(result)
  }
}








