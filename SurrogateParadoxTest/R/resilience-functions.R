################################################################################
####################### Functions to make the intervals ########################
################################################################################

gaussian_kernel <- function(x) {
  exp(-x^2 / 2) / sqrt(2 * pi)
}

smoother_fitter_extrapolate <- function(X, Y, kernel = gaussian_kernel, h) {
  # Precompute max Y at the left and right ends
  left_Y_max <- max(Y[X == min(X)])
  right_Y_max <- max(Y[X == max(X)])
  
  smoother <- function(x) {
    weights <- sapply(X, function(val) kernel((x - val)/h) / h)
    total_weight <- sum(weights)
    
    # Catch division by zero just in case
    if (total_weight == 0) {
      if (x < min(X)) {
        return(left_Y_max)
      } else if (x > max(X)) {
        return(right_Y_max)
      }
    }
    sum(weights * Y) / total_weight
  }
  smoother
}

gaussian_process_interval <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, 
                                    sigma2, theta, n.iter = 500,
                                    M = 100, q_quant = 0.1, plot = FALSE,
                                    intervals = TRUE, get_var = TRUE) {
  # Helper functions
  rbf_kernel <- function(x1, x2, length_scale = 1.0, variance = 1.0) {
    sqdist <- outer(x1, x2, function(a, b) (a - b)^2)
    variance * exp(-0.5 * sqdist / length_scale^2)
  }
  
  # First step is to estimate the mean function, which we will do using 
  # a kernel smoother
  mu0 <- smoother_fitter_extrapolate(
    s0.A, y0.A, h = bw.nrd(s0.A)
  )
  mu1 <- smoother_fitter_extrapolate(
    s1.A, y1.A, h = bw.nrd(s1.A)
  )
  
  # Mean parameter for MVN
  mean0 <- sapply(s0.B, mu0)
  mean1 <- sapply(s1.B, mu1)
  
  # Covariance parameter for MVN
  K0 <- rbf_kernel(s0.B, s0.B, theta, sigma2)
  K1 <- rbf_kernel(s1.B, s1.B, theta, sigma2)
  
  # Vectorized version of y_hat generation. Each row is one iter
  y0_hats <- mvrnorm(n.iter, mu = mean0, Sigma = K0)
  y1_hats <- mvrnorm(n.iter, mu = mean1, Sigma = K1)
  EY0s <- rowMeans(y0_hats)
  EY1s <- rowMeans(y1_hats)
  Delta_hats <- EY1s - EY0s
  
  if (get_var) {
    p_hats <- rep(NA, M)
    q_hats <- rep(NA, M)
    for (m in 1:M) {
      # Resample Study A
      control_idx <- sample(1:length(s0.A), length(s0.A), replace = TRUE)
      treat_idx <- sample(1:length(s1.A), length(s1.A), replace = TRUE)
      new_s0.A <- s0.A[control_idx]
      new_y0.A <- y0.A[control_idx]
      new_s1.A <- s1.A[treat_idx]
      new_y1.A <- y1.A[treat_idx]
      
      # New kernel smoothed Study A
      new_mu0 <- smoother_fitter_extrapolate(
        new_s0.A, new_y0.A, h = bw.nrd(s0.A)
      )
      new_mu1 <- smoother_fitter_extrapolate(
        new_s1.A, new_y1.A, h = bw.nrd(s1.A)
      )
      
      # Resample indexes for Study B
      control_idx_B <- sample(1:length(s0.B), length(s0.B), replace = TRUE)
      treat_idx_B <- sample(1:length(s1.B), length(s1.B), replace = TRUE)
      
      new_s0.B <- s0.B[control_idx_B]
      new_s1.B <- s1.B[treat_idx_B]
      
      new_mean0 <- sapply(new_s0.B, new_mu0)
      new_K0 <- rbf_kernel(new_s0.B, new_s0.B, theta, sigma2)
      
      new_mean1 <- sapply(new_s1.B, new_mu1)
      new_K1 <- rbf_kernel(new_s1.B, new_s1.B, theta, sigma2)
      
      Y0_mat <- mvrnorm(n.iter, new_mean0, new_K0)  # n.iter x d matrix
      Y1_mat <- mvrnorm(n.iter, new_mean1, new_K1)  # n.iter x d matrix
      
      # Compute means across rows
      mean_Y0 <- rowMeans(Y0_mat)
      mean_Y1 <- rowMeans(Y1_mat)
      
      # Compute the bootstrapped Delta_Bs
      Delta_Bs <- mean_Y1 - mean_Y0
      
      p_hats[m] <- mean(Delta_Bs < 0)
      q_hats[m] <- quantile(Delta_Bs, q_quant)[[1]]
    }
    p_se <- sd(p_hats)
    q_se <- sd(q_hats)
  }
  
  # Estimated p_hat and q_hat
  p_hat <- mean(Delta_hats < 0)
  q_hat <- quantile(Delta_hats, q_quant)[[1]]

  # Return the stuff
  result <- list(
    Delta_hats = Delta_hats,
    Delta_estimate = mean(Delta_hats),
    p_hat = p_hat,
    q_hat = q_hat
  )
  
  # Add in SE if requested
  if (get_var) {
    result$p_se = p_se
    result$q_se = q_se
  }
  
  if (plot) {
    # Scatterplot of points of A
    control_data <- data.frame(s = s0.A, y = y0.A)
    
    control_p <- ggplot() +
      geom_point(data = control_data, aes(x = s, y = y), color = "red", size = 0.5)
    
    treat_data <- data.frame(s = s1.A, y = y1.A)
    
    treat_p <- ggplot() +
      geom_point(data = treat_data, aes(x = s, y = y), color = "blue", size = 0.5)
    
    #To not overcrowd the plot, I'm just going to add the first 100
    for (i in 1:min(n.iter, 100)) {
      new_control <- data.frame(s = s0.B, y = y0_hats[i,])
      new_treat <- data.frame(s = s1.B, y = y1_hats[i,])
      control_p <- control_p +
        stat_smooth(data = new_control, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
      treat_p <- treat_p +
        stat_smooth(data = new_treat, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
    }
    
    control_p <- control_p +
      stat_smooth(data = control_data, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s0.B))

    treat_p <- treat_p +
      stat_smooth(data = treat_data, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s1.B))
    
    result$control_plot <- control_p
    result$treatment_plot <- treat_p
  }
  
  # Add intervals if requested
  if (get_var && intervals) {
    # Get the percentiles of the ps and the qs
    p_interval <- as.vector(quantile(p_hats, c(0.025, 0.975)))
    result$p_interval <- p_interval
    q_interval <- p_interval <- as.vector(quantile(q_hats, c(0.025, 0.975)))
    result$q_interval <- q_interval
  }
  return(result)
}

fourier_interval <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, var_vec, period,
                                     n.iter = 500, M = 100, q_quant = 0.1, plot = FALSE, 
                                     intervals = TRUE, get_var = TRUE) {
  
  mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
  mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))
  
  #predictions using mu functions from study A
  preds_0 <- sapply(s0.B, mu_hat_0)
  preds_1 <- sapply(s1.B, mu_hat_1)
  
  B0 <- 2 * pi / (period * (max(s0.A) - min(s0.A)))
  B1 <- 2 * pi / (period * (max(s1.A) - min(s1.A)))
  
  X0 <- cbind(preds_0, 
              1,
              sapply(s0.B, function(x) sin((x - min(s0.A)) / B0[1]) + cos((x - min(s0.A)) / B0[1])),
              sapply(s0.B, function(x) sin((x - min(s0.A)) / B0[2]) + cos((x - min(s0.A)) / B0[2])),
              sapply(s0.B, function(x) sin((x - min(s0.A)) / B0[3]) + cos((x - min(s0.A)) / B0[3]))
              )
  X1 <- cbind(preds_1, 
              1,
              sapply(s1.B, function(x) sin((x - min(s1.A)) / B1[1]) + cos((x - min(s1.A)) / B1[1])),
              sapply(s1.B, function(x) sin((x - min(s1.A)) / B1[2]) + cos((x - min(s1.A)) / B1[2])),
              sapply(s1.B, function(x) sin((x - min(s1.A)) / B1[3]) + cos((x - min(s1.A)) / B1[3]))
              )
  
  beta0matA <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
  beta1matA <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
  
  preds_0_varied <- beta0matA %*% t(X0)
  preds_1_varied <- beta1matA %*% t(X1)
  
  Deltas <- rowMeans(preds_1_varied) - rowMeans(preds_0_varied)

  q <- quantile(Deltas, q_quant)[[1]]
  
  # Now to do the part with the variance
  if (get_var) {
    p_hats <- rep(NA, M)
    q_hats <- rep(NA, M)
    for (m in 1:M) {
      # Resample Study A basically get a new function.
      control_idx <- sample(1:length(s0.A), length(s0.A), replace = TRUE)
      treat_idx <- sample(1:length(s1.A), length(s1.A), replace = TRUE)
      
      new_s0.A <- s0.A[control_idx]
      new_y0.A <- y0.A[control_idx]
      new_s1.A <- s1.A[treat_idx]
      new_y1.A <- y1.A[treat_idx]
      
      # New kernel smoothed Study A
      new_mu_hat_0 <- smoother_fitter_extrapolate(
        new_s0.A, new_y0.A, h = bw.nrd(s0.A)
      )
      new_mu_hat_1 <- smoother_fitter_extrapolate(
        new_s1.A, new_y1.A, h = bw.nrd(s1.A)
      )
      
      # Resample indexes for Study B
      control_idx_B <- sample(1:length(s0.B), length(s0.B), replace = TRUE)
      treat_idx_B <- sample(1:length(s1.B), length(s1.B), replace = TRUE)
      
      new_s0.B <- s0.B[control_idx_B]
      new_s1.B <- s1.B[treat_idx_B]
      
      # Now resample some Study B y values
      new_deltas <- rep(NA, n.iter)
      
      beta0matB <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
      beta1matB <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
      
      new_Y0.B <- beta0matB %*% t(cbind(sapply(new_s0.B, new_mu_hat_0), X0[,-1]))
      new_Y1.B <- beta1matB %*% t(cbind(sapply(new_s1.B, new_mu_hat_1), X1[,-1])) 
      
      new_deltas <- rowMeans(new_Y1.B) - rowMeans(new_Y0.B)

      p_hats[m] <- mean(new_deltas < 0)
      q_hats[m] <- quantile(new_deltas, q_quant)[[1]]
    }
    p_se <- sd(p_hats)
    q_se <- sd(q_hats)
  }
  
  # Return the stuff
  result <- list(
    Delta_hats = Deltas,
    Delta_estimate = mean(Deltas),
    p_hat = mean(Deltas < 0),
    q_hat = q
  )
  
  if (get_var) {
    result$p_se = p_se
    result$q_se = q_se
  }
  
  if (plot) {
    # Scatterplot of points of A
    control_data <- data.frame(s = s0.A, y = y0.A)
    
    control_p <- ggplot() +
      geom_point(data = control_data, aes(x = s, y = y), color = "red", size = 0.5)
    
    treat_data <- data.frame(s = s1.A, y = y1.A)
    
    treat_p <- ggplot() +
      geom_point(data = treat_data, aes(x = s, y = y), color = "blue", size = 0.5)
    
    #To not overcrowd the plot, I'm just going to add the first 100
    for (i in 1:min(n.iter, 100)) {
      new_control <- data.frame(s = s0.B, y = preds_0_varied[i,])
      new_treat <- data.frame(s = s1.B, y = preds_1_varied[i,])
      control_p <- control_p +
        stat_smooth(data = new_control, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
      treat_p <- treat_p +
        stat_smooth(data = new_treat, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
    }
    
    control_p <- control_p +
      stat_smooth(data = control_data, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s0.B))
    
    treat_p <- treat_p +
      stat_smooth(data = treat_data, aes(x = s, y = y), method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s1.B))
    
    result$control_plot <- control_p
    result$treatment_plot <- treat_p
  }
  return(result)
}

polynomial_interval <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, 
                                var_vec, n.iter = 500, M = 100, 
                                q_quant = 0.1, plot = FALSE, intervals = TRUE,
                                get_var = TRUE) {
  
  std.s0.B <- (s0.B - mean(s0.A)) / sd(s0.A)
  std.s1.B <- (s1.B - mean(s1.A)) / sd(s1.A)
  
  mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
  mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))
  
  mu1_result <- sapply(s1.B, mu_hat_1)
  mu0_result <- sapply(s0.B, mu_hat_0)
  
  X0 <- cbind(mu0_result, 1, std.s0.B, std.s0.B^2, std.s0.B^3)
  X1 <- cbind(mu1_result, 1, std.s1.B, std.s1.B^2, std.s1.B^3)
  
  beta0matA <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
  beta1matA <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))

  preds_0_varied <- beta0matA %*% t(X0)
  preds_1_varied <- beta1matA %*% t(X1)
  
  Deltas <- rowMeans(preds_1_varied) - rowMeans(preds_0_varied)
  
  # Now to do the part with the variance
  if (get_var) {
    p_hats <- rep(NA, M)
    q_hats <- rep(NA, M)
    for (m in 1:M) {
      # Resample Study A basically get a new function.
      control_idx <- sample(1:length(s0.A), length(s0.A), replace = TRUE)
      treat_idx <- sample(1:length(s1.A), length(s1.A), replace = TRUE)

      new_s0.A <- s0.A[control_idx]
      new_y0.A <- y0.A[control_idx]
      
      new_s1.A <- s1.A[treat_idx]
      new_y1.A <- y1.A[treat_idx]
      
      # New kernel smoothed Study A
      new_mu_hat_0 <- smoother_fitter_extrapolate(
        new_s0.A, new_y0.A, h = bw.nrd(s0.A)
      )
      new_mu_hat_1 <- smoother_fitter_extrapolate(
        new_s1.A, new_y1.A, h = bw.nrd(s1.A)
      )
      
      # Resample indexes for Study B
      control_idx_B <- sample(1:length(s0.B), length(s0.B), replace = TRUE)
      treat_idx_B <- sample(1:length(s1.B), length(s1.B), replace = TRUE)
      
      new_s0.B <- s0.B[control_idx_B]
      new_std.s0.B <- (new_s0.B - mean(s0.A)) / sd(s0.A)
      new_s1.B <- s1.B[treat_idx_B]
      new_std.s1.B <- (new_s1.B - mean(s1.A)) / sd(s1.A)
      
      # Now resample some Study B y values
      beta0matB <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
      beta1matB <- cbind(1, mvrnorm(n = n.iter, mu = rep(0, length(var_vec)), Sigma = diag(var_vec)))
      
      new_mu0_base <- sapply(new_s0.B, new_mu_hat_0)
      new_mu1_base <- sapply(new_s1.B, new_mu_hat_1)
      
      new_X0 <- cbind(new_mu0_base, 1, new_std.s0.B, new_std.s0.B^2, new_std.s0.B^3)
      new_X1 <- cbind(new_mu1_base, 1, new_std.s1.B, new_std.s1.B^2, new_std.s1.B^3)
      
      new_y0.B <- beta0matB %*% t(new_X0)
      new_y1.B <- beta1matB %*% t(new_X1)
      
      new_deltas <- rowMeans(new_y1.B) - rowMeans(new_y0.B)
      p_hats[m] <- mean(new_deltas < 0)
      q_hats[m] <- quantile(new_deltas, q_quant)[[1]]
    }
    p_se <- sd(p_hats)
    q_se <- sd(q_hats)
  }
  
  # Get q
  q <- quantile(Deltas, q_quant)[[1]]
  
  # Return the stuff
  result <- list(
    Delta_hats = Deltas,
    Delta_estimate = mean(Deltas),
    p_hat = mean(Deltas < 0),
    q_hat = q
  )
  
  if (plot) {
    # Scatterplot of points of A
    control_data <- data.frame(s = s0.A, y = y0.A)
    
    control_p <- ggplot() +
      geom_point(data = control_data, aes(x = s, y = y), color = "red", size = 0.5)
    
    treat_data <- data.frame(s = s1.A, y = y1.A)
    
    treat_p <- ggplot() +
      geom_point(data = treat_data, aes(x = s, y = y), color = "blue", size = 0.5)
    
    #To not overcrowd the plot, I'm just going to add the first 100
    for (i in 1:min(n.iter, 100)) {
      new_control <- data.frame(s = s0.B, y = preds_0_varied[i,])
      new_treat <- data.frame(s = s1.B, y = preds_1_varied[i,])
      control_p <- control_p +
        stat_smooth(data = new_control, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
      treat_p <- treat_p +
        stat_smooth(data = new_treat, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "gray", size = 0.5, alpha = 0.2, se = FALSE)
    }
    
    control_p <- control_p +
      stat_smooth(data = control_data, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s0.B))
    
    treat_p <- treat_p +
      stat_smooth(data = treat_data, aes(x = s, y = y),method = "loess", formula = y ~ x, color = "black", se = FALSE) +
      coord_cartesian(xlim = range(s1.B))
    
    result$control_plot <- control_p
    result$treatment_plot <- treat_p
  }
  
  if (get_var) {
    result$p_se = p_se
    result$q_se = q_se
  }
  
  return(result)
}

################################################################################
########################## Functions to generate data ##########################
################################################################################

rbf_kernel <- function(x1, x2, theta, sigma2) {
  # Compute the kernel matrix
  sqdist <- outer(x1, x2, function(a, b) (a - b)^2)
  sigma2 * exp(-0.5 * sqdist / theta^2)
}


################################################################################
######################## Functions for resilience set ##########################
################################################################################

gp_resilience_set <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, 
                              sigma2_vals, theta_vals, alpha = 0.05) {
  
  mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
  mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))

  mu0 <- sapply(s0.B, mu_hat_0)
  mu1 <- sapply(s1.B, mu_hat_1)
  n0 <- length(mu0)
  n1 <- length(mu1)
  a0 <- rep(1 / n0, n0)
  a1 <- rep(1 / n1, n1)
  
  f <- function(theta, sigma2) {
    K0 <- rbf_kernel(s0.B, s0.B, theta, sigma2)
    K1 <- rbf_kernel(s1.B, s1.B, theta, sigma2)
    
    DeltaB_mean <- t(a1) %*% mu1 - t(a0) %*% mu0
    DeltaB_variance <- t(a1) %*% K1 %*% a1 + t(a0) %*% K0 %*% a0
    
    return(pnorm(0, mean = DeltaB_mean, sd = sqrt(DeltaB_variance)) - alpha)
  }
  
  # Plot the feasible region and the f(x) = 0 contour line using a coarse grid just for visualization
  plot_grid <- expand.grid(
    x1 = theta_vals,
    x2 = sigma2_vals
  )
  plot_grid$z <- mapply(f, plot_grid$x1, plot_grid$x2)
  
  plot <- ggplot() +
    geom_tile(data = subset(plot_grid, z <= 0),
              aes(x = x1, y = x2), fill = "lightblue", alpha = 0.5) +
    geom_contour(data = plot_grid, aes(x = x1, y = x2, z = z),
                 breaks = 0, color = "black", linewidth = 1) +
    labs(
      x = expression(theta),
      y = expression(sigma^2)
    ) +
    theme_minimal()
  return(plot)
}

polynomial_resilience_set <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B,
                                      sig1_values, sig2_values, alpha = 0.05) {
  
  mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
  mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))
  
  # Apply values 
  mu0 <- sapply(s0.B, mu_hat_0)
  mu1 <- sapply(s1.B, mu_hat_1)
  n0 <- length(mu0)
  n1 <- length(mu1)
  a0 <- rep(1 / n0, n0)
  a1 <- rep(1 / n1, n1)
  # To make the M matrix
  f0_1 <- function(s) 1
  f0_2 <- function(s) (s - mean(s0.A)) / sd(s0.A)
  f0_3 <- function(s) ((s - mean(s0.A)) / sd(s0.A))^2
  f0_4 <- function(s) ((s - mean(s0.A)) / sd(s0.A))^3
  funcs_0 <- list(f0_1, f0_2, f0_3, f0_4)
  f1_1 <- function(s) 1
  f1_2 <- function(s) (s - mean(s1.A)) / sd(s1.A)
  f1_3 <- function(s) ((s - mean(s1.A)) / sd(s1.A))^2
  f1_4 <- function(s) ((s - mean(s1.A)) / sd(s1.A))^3
  funcs_1 <- list(f1_1, f1_2, f1_3, f1_4)
  
  M0 <- t(sapply(s0.B, function(x) sapply(funcs_0, function(f) f(x))))
  M1 <- t(sapply(s1.B, function(x) sapply(funcs_1, function(f) f(x))))
  
  f_poly <- function(x, y) {
    Sigma <- diag(c(x, x, y, y))
    DeltaB_mean <- t(a1) %*% mu1 - t(a0) %*% mu0
    DeltaB_variance <- t(a1) %*% M1 %*% Sigma %*% t(M1) %*% a1 + 
      t(a0) %*% M0 %*% Sigma %*% t(M0) %*% a0
    return(pnorm(0, mean = DeltaB_mean, sd = sqrt(DeltaB_variance)) - alpha)
  }

  # Plot the feasible region and the f(x) = 0 contour line using a coarse grid just for visualization
  plot_grid <- expand.grid(x = sig1_values, y = sig2_values)
  plot_grid$z <- mapply(f_poly, plot_grid$x, plot_grid$y)
  plot_grid$z_bin <- ifelse(plot_grid$z < 0, "lightblue", "white")
  
  plot <- ggplot() +
    geom_tile(data = plot_grid,
              aes(x = x, y = y, fill = z_bin), alpha = 0.5) +
    scale_fill_identity() +
    geom_contour(data = plot_grid, aes(x = x, y = y, z = z),
                 breaks = 0, color = "black", linewidth = 1) +
    labs(
      x = expression(sigma[11]^2),
      y = expression(sigma[22]^2)
    ) +
    theme_minimal()
  return(plot)
}

fourier_resilience_set <- function(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, 
                                   sig1_values, sig2_values, alpha = 0.05) {
  
  mu_hat_1 <- smoother_fitter_extrapolate(s1.A, y1.A, h = bw.nrd(s1.A))
  mu_hat_0 <- smoother_fitter_extrapolate(s0.A, y0.A, h = bw.nrd(s0.A))
  
  # Apply values 
  mu0 <- sapply(s0.B, mu_hat_0)
  mu1 <- sapply(s1.B, mu_hat_1)
  n0 <- length(mu0)
  n1 <- length(mu1)
  a0 <- rep(1 / n0, n0)
  a1 <- rep(1 / n1, n1)
  
  # Some necessary parameters
  period = c(0.5, 0.25, 0.1, 0.05)
  B0 <- 2 * pi / (period * (max(y0.A) - min(y0.A)))
  B1 <- 2 * pi / (period * (max(y1.A) - min(y1.A)))
  
  # To make the M matrix
  f0_1 <- function(s) 1
  f0_2 <- function(s) sin((s - min(y0.A)) / B0[1]) + cos((s - min(y0.A)) / B0[1])
  f0_3 <- function(s) sin((s - min(y0.A)) / B0[2]) + cos((s - min(y0.A)) / B0[2])
  f0_4 <- function(s) sin((s - min(y0.A)) / B0[3]) + cos((s - min(y0.A)) / B0[3])
  funcs_0 <- list(f0_1, f0_2, f0_3, f0_4)
  f1_1 <- function(s) 1
  f1_2 <- function(s) sin((s - min(y1.A)) / B1[1]) + cos((s - min(y1.A)) / B1[1])
  f1_3 <- function(s) sin((s - min(y1.A)) / B1[2]) + cos((s - min(y1.A)) / B1[2])
  f1_4 <- function(s) sin((s - min(y1.A)) / B1[3]) + cos((s - min(y1.A)) / B1[3])
  funcs_1 <- list(f1_1, f1_2, f1_3, f1_4)
  
  M0 <- t(sapply(s0.B, function(x) sapply(funcs_0, function(f) f(x))))
  M1 <- t(sapply(s1.B, function(x) sapply(funcs_1, function(f) f(x))))
  
  f_fourier <- function(x, y) {
    Sigma <- diag(c(x, x, y, y))
    DeltaB_mean <- t(a1) %*% mu1 - t(a0) %*% mu0
    DeltaB_variance <- t(a1) %*% M1 %*% Sigma %*% t(M1) %*% a1 + 
      t(a0) %*% M0 %*% Sigma %*% t(M0) %*% a0
    return(pnorm(0, mean = DeltaB_mean, sd = sqrt(DeltaB_variance)) - alpha)
  }
  
  # Plot the feasible region and the f(x) = 0 contour line using a coarse grid just for visualization
  plot_grid <- expand.grid(
    x1 = sig1_values,
    x2 = sig2_values
  )
  plot_grid$z <- mapply(f_fourier, plot_grid$x1, plot_grid$x2)
  plot_grid$z_bin <- ifelse(plot_grid$z < 0, "lightblue", "white")
  
  ggplot() +
    geom_tile(data = plot_grid,
              aes(x = x1, y = x2, fill = z_bin), alpha = 0.5) +
    scale_fill_identity() +
    geom_contour(data = plot_grid, aes(x = x1, y = x2, z = z),
                 breaks = 0, color = "black", linewidth = 1) +
    labs(
      x = expression(sigma[11]^2),
      y = expression(sigma[22]^2)
    ) +
    theme_minimal()
}

