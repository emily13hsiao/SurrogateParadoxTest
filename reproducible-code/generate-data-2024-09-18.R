# Generate data for settings in paper.

###############################################################################
###########################      Assumption 1      ############################
###############################################################################

a1.data <- function(n_sample, setting) {
  if (setting == 1) {
    s1 <- rnorm(n_sample, mean = 24, sd = 1)
    s0 <- rnorm(n_sample, mean = 23, sd = 1)
  } else if (setting == 2) {
    s1 <- rexp(n_sample)
    s0 <- rgamma(n_sample, 0.25)
  } else if (setting == 3) {
    s1 <- rnorm(n_sample, mean = 10, sd = 1)
    s0 <- rnorm(n_sample, mean = 1, sd = 1)
  } else if (setting == 4) {
    s1 <- rnorm(n_sample)
    s0 <- rnorm(n_sample)
  } else if (setting == 5) {
    s1 <- rnorm(n_sample, mean = 0, sd = sqrt(5))
    s0 <- rnorm(n_sample, mean = -2, sd = sqrt(0.5))
  } else if (setting == 6) {
    s1 <- rnorm(n_sample, mean = 0, sd = sqrt(5))
    s0 <- rnorm(n_sample, mean = 0, sd = 1)
  } else if (setting == 7) {
    s1 <- rnorm(n_sample, mean = 0, sd = sqrt(1.5))
    s0 <- rnorm(n_sample, mean = 0, sd = 1)
  } else if (setting == 8) {
    s1 <- rnorm(n_sample, mean = 0, sd = sqrt(5))
    s0 <- rnorm(n_sample, mean = 1, sd = sqrt(2))
  } else if (setting == 9) {
    s1 <- rnorm(n_sample, mean = 0, sd = 1)
    s0 <- rnorm(n_sample, mean = 1, sd = 1)
  }
  return(list(s0 = s0, s1 = s1))
}

###############################################################################
###########################      Assumption 2      ############################
###############################################################################

g_1 <- function(x, M) {
  if (x <= 0.5) {
    15 * (x - 0.5)^3 + M * (x - 0.5)
  } else {
    M * (x - 0.5)
  }
}

g_2 <- function(x) exp(-250 * (x - 0.25)^2)

make_g <- function(M) {
  f <- function(x) g_1(x, M) - g_2(x)
  return(f)
}

a2.data = function(n_sample, setting, sim_sd) {
  if (setting == "decreasing") {
    x <- runif(n_sample)
    y=-2*x + rnorm(length(x), 0, sim_sd)
  }
  if(setting == "flat") {
    x <- runif(n_sample)
    y=rnorm(length(x), 0, sim_sd)
  }
  if(setting == "increasing") {
    x <- runif(n_sample)
    y=x^2 + rnorm(length(x),0, sim_sd)
  }
  if(setting == "hall") {
    g <- make_g(0.3)
    x <- runif(n_sample)
    y <- sapply(x, function(x_i) rnorm(1, g(x_i), sim_sd))
  }
  if (setting == "parabola") {
    x <- runif(n_sample)
    y = (x - 0.5)^2 + rnorm(length(x), 0, sim_sd)
  } 
  if (setting == "increasing_linear") {
    x <- runif(n_sample)
    y <- 2*x + rnorm(length(x), 0, sim_sd)
  }
  return(list("x" = x, "y" = y))
}

###############################################################################
###########################      Assumption 3      ############################
###############################################################################

a3.data <- function(n_sample, setting) {
  if (setting == 1) {
    m_c <- function(s) s
    m_t <- function(s) 3 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 2) {
    m_c <- function(s) s
    m_t <- function(s) 3 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 3) {
    m_c <- function(s) 1 + 2 * s
    m_t <- function(s) 1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 4) {
    m_c <- function(s) 1 + 2 * s
    m_t <- function(s) 1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 5) {
    m_c <- function(s) s
    m_t <- function(s) 1.5 + s^2
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 3))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 3))
  } else if (setting == 6) {
    m_c <- function(s) 3 + 2 * s / 3
    m_t <- function(s) -1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 7) {
    m_c <- function(s) 3 + 2 * s / 3
    m_t <- function(s) -1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 8) {
    m_c <- function(s) 2 * s 
    m_t <- function(s) s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 9) {
    m_c <- function(s) s 
    m_t <- function(s) 0.5 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  }
  return(list("s0" = s_c, "y0" = y_c, "s1" = s_t, "y1" = y_t))
}

generate.data <- function(n_sample, setting) {
  if (setting == 1) {
    m_c <- function(s) s
    m_t <- function(s) 3 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 2) {
    m_c <- function(s) s
    m_t <- function(s) 3 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 3) {
    m_c <- function(s) 1 + 2 * s
    m_t <- function(s) 1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 4) {
    m_c <- function(s) 1 + 2 * s
    m_t <- function(s) 1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 5) {
    m_c <- function(s) s
    m_t <- function(s) 1.5 + s^2
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 3))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 3))
  } else if (setting == 6) {
    m_c <- function(s) 3 + 2 * s / 3
    m_t <- function(s) -1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 7) {
    m_c <- function(s) 3 + 2 * s / 3
    m_t <- function(s) -1 + 2 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 5))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 5))
  } else if (setting == 8) {
    m_c <- function(s) 2 * s 
    m_t <- function(s) s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  } else if (setting == 9) {
    m_c <- function(s) s 
    m_t <- function(s) 0.5 * s
    
    s_c <- rnorm(n_sample, 3, 1)
    y_c <- sapply(s_c, function(s) rnorm(1, m_c(s), 1))
    s_t <- rnorm(n_sample, 3, 1)
    y_t <- sapply(s_t, function(s) rnorm(1, m_t(s), 1))
  }
  return(list("s0" = s_c, "y0" = y_c, "s1" = s_t, "y1" = y_t))
}

