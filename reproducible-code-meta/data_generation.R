# Data generating code
library(MASS)
source("functions-with-se.R")

setting1 <- function(n.studies, n.k) {
  # This is linear setting, high prob, with overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 1
  v2 <- 1
  
  # True mean functions
  m0 <- function(s) 2 * s - 1
  m1 <- function(s) s + 3
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A <- sort(rnorm(n.k, 3, sqrt(3)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A <- sort(rnorm(n.k, 4, sqrt(3)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, 4.75, sqrt(1)))
  s1.B <- sort(rnorm(n.k, 5.25, sqrt(1)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}

setting2 <- function(n.studies, n.k) {
  # This is polynomial setting, high prob, with overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 1
  v2 <- 1
  
  # True mean functions
  m0 <- function(s) (s - 0.5)^2 - 1
  m1 <- function(s) 3 * s + 1
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A <- sort(rnorm(n.k, 0.9, sqrt(1.5)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A <- sort(rnorm(n.k, 2.2, sqrt(4.5)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, -0.7, sqrt(1)))
  s1.B <- sort(rnorm(n.k, -0.2, sqrt(2)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}

setting3 <- function(n.studies, n.k) {
  # This is Fourier setting, high prob, overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 1
  v2 <- 1
  
  # True mean functions
  m0 <- function(s) 0.2 + 0.4 * sin(s) + 0.4 * cos(s)
  m1 <- function(s) 0.6 + 0.85 * sin(s) + 0.85 * cos(s)
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A<-sort(rnorm(n.k, 5, sqrt(1)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A<-sort(rnorm(n.k, 6, sqrt(2)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, 4.1, sqrt(0.5)))
  s1.B <- sort(rnorm(n.k, 4.5, sqrt(0.5)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}


setting4 <- function(n.studies, n.k) {
  # This is low probability, linear setting, overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 1
  v2 <- 1
  
  # True mean functions
  m0 <- function(s) 1.5 * s + 1 
  m1 <- function(s) 3 * s - 2
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A <- sort(rnorm(n.k, 2, sqrt(3)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A <- sort(rnorm(n.k, 3, sqrt(3)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, 1.75, sqrt(1)))
  s1.B <- sort(rnorm(n.k, 2.75, sqrt(1)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}

setting5 <- function(n.studies, n.k) {
  # This is polynomial setting, low prob, with overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 1
  v2 <- 1
  
  # True mean functions
  m0 <- function(s) (s - 0.5)^2-1 
  m1 <- function(s) 3 * s + 1
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A <- sort(rnorm(n.k, 0.9, sqrt(1.5)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A <- sort(rnorm(n.k, 2.2, sqrt(4.5)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, -0.08, sqrt(1)))
  s1.B <- sort(rnorm(n.k, 0.45, sqrt(2)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}

setting6 <- function(n.studies, n.k) {
  # This is Fourier setting, low prob, with overlapping support
  
  # Other parameters
  theta <- 5
  sig2 <- 0.1
  v2 <- 0.5
  
  # True mean functions
  m0 <- function(s) 0.2 + 0.4 * sin(s) + 0.4 * cos(s)
  m1 <- function(s) 0.6 + 0.85 * sin(s) + 0.85 * cos(s)
  
  # Generate more studies
  studies <- list()
  for (i in 1:n.studies) {
    # Study A Data
    s0.A <- sort(rnorm(n.k, 5, sqrt(1)))
    y0.A <- mvrnorm(1, mu = sapply(s0.A, m0), Sigma = rbf_kernel(s0.A, s0.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    s1.A <- sort(rnorm(n.k, 6, sqrt(2)))
    y1.A <- mvrnorm(1, mu = sapply(s1.A, m1), Sigma = rbf_kernel(s1.A, s1.A, theta, sig2)) +
      rnorm(n.k, 0, sqrt(v2))
    
    studies[[i]] <- list(s0.A = s0.A, y0.A = y0.A, s1.A = s1.A, y1.A = y1.A)
  }
  s0.B <- sort(rnorm(n.k, 5.5, sqrt(0.5)))
  s1.B <- sort(rnorm(n.k, 6.5, sqrt(0.5)))
  y0.B <- mvrnorm(1, mu = sapply(s0.B, m0), Sigma = rbf_kernel(s0.B, s0.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  y1.B <- mvrnorm(1, mu = sapply(s1.B, m1), Sigma = rbf_kernel(s1.B, s1.B, theta, sig2)) +
    rnorm(n.k, 0, sqrt(v2))
  return(list(data = studies, s0.B = s0.B, s1.B = s1.B, y0.B = y0.B, y1.B = y1.B))
}

#generate data function 

generate_data <- function(setting, n.studies, n.k) {
  if (setting == 1) {
    data <- setting1(n.studies, n.k)
  } else if (setting == 2) {
    data <- setting2(n.studies, n.k)
  } else if (setting == 3) {
    data <- setting3(n.studies, n.k)
  } else if (setting == 4) {
    data <- setting4(n.studies, n.k)
  } else if (setting == 5) {
    data <- setting5(n.studies, n.k)
  } else if (setting == 6) {
    data <- setting6(n.studies, n.k)
  } 
  
  #flatten the data, take it out of lists
  df_list <- list()
  
  # Loop over the K studies
  for (i in seq_along(data$data)) {
    study <- data$data[[i]]
    
    # control arm (G = 0)
    df0 <- data.frame(
      S = study$s0.A,
      Y = study$y0.A,
      G = 0,
      study = i
    )
    
    # treated arm (G = 1)
    df1 <- data.frame(
      S = study$s1.A,
      Y = study$y1.A,
      G = 1,
      study = i
    )
    
    df_list[[length(df_list) + 1]] <- df0
    df_list[[length(df_list) + 1]] <- df1
  }
  
  
  # combine everything
  flat_df <- do.call(rbind, df_list)
  rownames(flat_df) <- NULL
  
  return(list(data = flat_df, s0.B = data$s0.B, s1.B = data$s1.B, y0.B = data$y0.B, y1.B = data$y1.B))
}