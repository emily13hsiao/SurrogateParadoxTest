
# Functions with the additional bootstrap

library(tidyverse)
library(MASS)
library(splines)
library(Rsurrogate)
library(RcppArmadillo)
library(Rcpp)
library(mvmeta)
library(parallel)

#################################################
## Functions for meta-analytic resilience paper
#################################################

#RBF kernel, uses C++

Rcpp::cppFunction(depends = "RcppArmadillo", code = '
arma::mat rbf_kernel(const arma::vec& x1, const arma::vec& x2, double length_scale, double variance) {
  int n1 = x1.n_elem;
  int n2 = x2.n_elem;
  arma::mat K(n1, n2);
  
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      double diff = x1[i] - x2[j];
      K(i, j) = variance * exp(-0.5 * (diff * diff) / (length_scale * length_scale));
    }
  }
  
  return K;
}
')

# Matrix inversion and likelihood contribution, C++

cppFunction(depends = "RcppArmadillo", code = '
double loglik_contrib(const arma::mat& Ck, const arma::vec& yk, const arma::vec& mk) {
  arma::vec v = yk - mk;
  arma::vec sol = arma::solve(Ck, v);
  double quad = arma::as_scalar(v.t() * sol);
  double logdet_val = 0.0;
  double sign = 0.0;
  arma::log_det(logdet_val, sign, Ck);
  return 0.5 * (logdet_val + quad);
}
')

# Just calculates the probability
calculate_p_hat0 <- function(data, s0.B, s1.B,
                          degree = 3,
                          n.reps = 200,
                          spline = FALSE) {
  
  # set up as list
  unique.studies <- unique(data$study)
  
  studies <- vector("list", length(unique.studies))
  for (j in seq_along(unique.studies)) {
    s0.A <- data$S[data$study == unique.studies[j] & data$G == 0]
    s1.A <- data$S[data$study == unique.studies[j] & data$G == 1]
    y0.A <- data$Y[data$study == unique.studies[j] & data$G == 0]
    y1.A <- data$Y[data$study == unique.studies[j] & data$G == 1]
    studies[[j]] <- list(
      i    = unique.studies[j],
      s0.A = s0.A, y0.A = y0.A,
      s1.A = s1.A, y1.A = y1.A
    )
  }
  
  ## ------------------------------------------------------------------------
  ## Likelihood definitions 
  ## ------------------------------------------------------------------------
  
  if (!spline) {
    # negative log likelihood without spline
    negloglike <- function(par, studies, treat) {
      sigma2 <- par[1]
      theta  <- par[2]
      v2     <- par[3]
      beta   <- par[-c(1:3)]
      
      loglik <- 0
      for (k in seq_along(studies)) {
        if (treat) {
          xk <- studies[[k]]$s1.A
          yk <- studies[[k]]$y1.A
        } else {
          xk <- studies[[k]]$s0.A
          yk <- studies[[k]]$y0.A
        }
        
        xbasis <- poly(xk, degree = degree, raw = TRUE)
        mk     <- cbind(1, xbasis) %*% beta
        
        Ck <- rbf_kernel(xk, xk, theta, sigma2) + v2 * diag(length(xk))
        loglik <- loglik + loglik_contrib(Ck, yk, mk)
      }
      loglik
    }
    
  } else {
    control.vec <- c(s0.B, data$S[data$G == 0])
    treat.vec   <- c(s1.B, data$S[data$G == 1])
    boundary_knots0 <- c(min(control.vec) - 1, max(control.vec) + 1)
    boundary_knots1 <- c(min(treat.vec)   - 1, max(treat.vec)   + 1)
    knots0 <- as.vector(quantile(control.vec, c(0.33, 0.67))) 
    knots1 <- as.vector(quantile(treat.vec,   c(0.33, 0.67))) 
    
    negloglikespline <- function(par, studies, treat, knots, boundary_knots) {
      sigma2 <- par[1]
      theta  <- par[2]
      v2     <- par[3]
      beta   <- par[-c(1:3)]
      
      loglik <- 0
      for (k in seq_along(studies)) {
        if (treat) {
          xk <- studies[[k]]$s1.A
          yk <- studies[[k]]$y1.A
        } else {
          xk <- studies[[k]]$s0.A
          yk <- studies[[k]]$y0.A
        }
        xbasis <- bs(xk, knots = knots,
                     Boundary.knots = boundary_knots,
                     intercept = TRUE)
        mk <- xbasis %*% beta
        
        Ck <- rbf_kernel(xk, xk, theta, sigma2) + v2 * diag(length(xk))
        loglik <- loglik + loglik_contrib(Ck, yk, mk)
      }
      loglik
    }
  }
  
  ## ------------------------------------------------------------------------
  ## Original fit: estimate parameters, simulate Y, compute p
  ## ------------------------------------------------------------------------
  
  if (!spline) {
    initial_params <- rep(1, 3 + 1 + degree)
    
    opt0 <- nlminb(start = initial_params,
                   objective = negloglike,
                   studies   = studies,
                   treat     = FALSE,
                   lower     = c(c(0.1, 0.1, 0.1), rep(-Inf, 1 + degree)),
                   upper     = rep(Inf, 3 + 1 + degree))
    par0 <- opt0$par
    
    opt1 <- nlminb(start = initial_params,
                   objective = negloglike,
                   studies   = studies,
                   treat     = TRUE,
                   lower     = c(c(0.1, 0.1, 0.1), rep(-Inf, 1 + degree)),
                   upper     = rep(Inf, 3 + 1 + degree))
    par1 <- opt1$par
    
    sigma2_hat0 <- par0[1]
    theta_hat0  <- par0[2]
    v2_hat0     <- par0[3]
    beta0_hat   <- par0[-c(1:3)]
    mu0_hat     <- cbind(1, poly(s0.B, degree = degree, raw = TRUE)) %*% beta0_hat 
    
    sigma2_hat1 <- par1[1]
    theta_hat1  <- par1[2]
    v2_hat1     <- par1[3]
    beta1_hat   <- par1[-c(1:3)]
    mu1_hat     <- cbind(1, poly(s1.B, degree = degree, raw = TRUE)) %*% beta1_hat 
    
  } else {
    initial_params <- rep(1, length(knots0) + degree + 4)
    
    opt0 <- nlminb(start = initial_params,
                   objective = negloglikespline,
                   studies   = studies,
                   treat     = FALSE,
                   knots     = knots0,
                   boundary_knots = boundary_knots0,
                   lower     = c(c(0.1, 0.1, 0.1),
                                 rep(-Inf, length(knots0) + degree + 1)),
                   upper     = rep(Inf, length(knots0) + degree + 4))
    par0 <- opt0$par
    
    opt1 <- nlminb(start = initial_params,
                   objective = negloglikespline,
                   studies   = studies,
                   treat     = TRUE,
                   knots     = knots1,
                   boundary_knots = boundary_knots1,
                   lower     = c(c(0.1, 0.1, 0.1),
                                 rep(-Inf, length(knots1) + degree + 1)),
                   upper     = rep(Inf, length(knots1) + degree + 4))
    par1 <- opt1$par
    
    sigma2_hat0 <- par0[1]
    theta_hat0  <- par0[2]
    v2_hat0     <- par0[3]
    beta0_hat   <- par0[-c(1:3)]
    mu0_hat     <- bs(s0.B, knots = knots0,
                      Boundary.knots = boundary_knots0,
                      intercept = TRUE) %*% beta0_hat
    
    sigma2_hat1 <- par1[1]
    theta_hat1  <- par1[2]
    v2_hat1     <- par1[3]
    beta1_hat   <- par1[-c(1:3)]
    mu1_hat     <- bs(s1.B, knots = knots1,
                      Boundary.knots = boundary_knots1,
                      intercept = TRUE) %*% beta1_hat
  }

  # Simulate Y given B-surrogate and fitted parameters
  Sigma0 <- rbf_kernel(s0.B, s0.B, theta_hat0, sigma2_hat0) +
    diag(v2_hat0, length(s0.B))
  Sigma1 <- rbf_kernel(s1.B, s1.B, theta_hat1, sigma2_hat1) +
    diag(v2_hat1, length(s1.B))
  
  Y0 <- mvrnorm(n.reps, mu = mu0_hat, Sigma = Sigma0)
  Y1 <- mvrnorm(n.reps, mu = mu1_hat, Sigma = Sigma1)
  
  Delta_Bs <- rowMeans(Y1) - rowMeans(Y0)
  p <- mean(Delta_Bs < 0)
  
  # Base return object (as you had)
  res <- list(
    p           = p,
    theta_hat0  = theta_hat0,
    theta_hat1  = theta_hat1,
    sigma2_hat0 = sigma2_hat0,
    sigma2_hat1 = sigma2_hat1,
    v2_hat0     = v2_hat0,
    v2_hat1     = v2_hat1,
    beta0_hat   = beta0_hat,
    beta1_hat   = beta1_hat,
    par0        = par0,
    par1        = par1
  )

  return(res)
}

# -------------------------------------------------------------------
# Helper function (same as before)
# -------------------------------------------------------------------
make_basis <- function(x, degree = 3, use_spline = FALSE, knots = NULL, boundary_knots = NULL) {
  if (!use_spline) {
    # Simple polynomial of given degree
    cbind(1, poly(x, degree = degree, raw = TRUE))
  } else {
    # Spline with specified knots and boundary knots
    if (is.null(knots)) stop("Knots must be provided for splines")
    if (is.null(boundary_knots)) stop("Boundary knots must be provided for splines")
    splines::bs(x, degree = degree, knots = knots,
                Boundary.knots = boundary_knots, intercept = TRUE)
  }
}

# -------------------------------------------------------------------
# Updated calculate_p_hat using make_basis()
# -------------------------------------------------------------------
calculate_p_hat <- function(data, s0.B, s1.B,
                            degree = 3,
                            n.reps = 200,
                            use_spline = FALSE, 
                            knots0, knots1,
                            boundary_knots0, boundary_knots1) {
  
  # set up study list
  unique.studies <- unique(data$study)
  studies <- lapply(unique.studies, function(st) {
    list(
      i    = st,
      s0.A = data$S[data$study == st & data$G == 0],
      y0.A = data$Y[data$study == st & data$G == 0],
      s1.A = data$S[data$study == st & data$G == 1],
      y1.A = data$Y[data$study == st & data$G == 1]
    )
  })
  
  # -------------------------------------------------------------------
  # Negative log-likelihood
  # -------------------------------------------------------------------
  negloglike <- function(par, studies, treat) {
    sigma2 <- par[1]
    theta  <- par[2]
    v2     <- par[3]
    beta   <- par[-c(1:3)]
    
    loglik <- 0
    for (k in seq_along(studies)) {
      if (treat) {
        xk <- studies[[k]]$s1.A
        yk <- studies[[k]]$y1.A
      } else {
        xk <- studies[[k]]$s0.A
        yk <- studies[[k]]$y0.A
      }
      
      Xk <- make_basis(xk, degree = degree, use_spline = use_spline,
                       knots = if(treat) knots1 else knots0,
                       boundary_knots = if(treat) boundary_knots1 else boundary_knots0)
      mk <- Xk %*% beta
      
      Ck <- rbf_kernel(xk, xk, theta, sigma2) + v2 * diag(length(xk))
      loglik <- loglik + loglik_contrib(Ck, yk, mk)
    }
    loglik
  }
  
  # -------------------------------------------------------------------
  # Fit parameters
  # -------------------------------------------------------------------
  p_dim <- if(use_spline) length(knots0) + degree + 1 else degree + 1
  initial_params <- rep(1, 3 + p_dim)
  
  opt0 <- nlminb(start = initial_params,
                 objective = negloglike,
                 studies   = studies,
                 treat     = FALSE,
                 lower     = c(c(0.1, 0.1, 0.1), rep(-Inf, p_dim)),
                 upper     = rep(Inf, 3 + p_dim))
  
  opt1 <- nlminb(start = initial_params,
                 objective = negloglike,
                 studies   = studies,
                 treat     = TRUE,
                 lower     = c(c(0.1, 0.1, 0.1), rep(-Inf, p_dim)),
                 upper     = rep(Inf, 3 + p_dim))
  
  par0 <- opt0$par
  par1 <- opt1$par
  beta0 <- par0[-c(1:3)]
  beta1 <- par1[-c(1:3)]
  
  # ---------------------------
  # REPLACE: use make_basis for mu
  # ---------------------------
  
  X0B <- make_basis(s0.B, degree = degree, use_spline = use_spline, knots = knots0, boundary_knots = boundary_knots0)
  X1B <- make_basis(s1.B, degree = degree, use_spline = use_spline, knots = knots1, boundary_knots = boundary_knots1)
  
  mu0_hat <- X0B %*% beta0
  mu1_hat <- X1B %*% beta1
  
  # -------------------------------------------------------------------
  # Simulate Y given B-surrogate and fitted parameters
  # -------------------------------------------------------------------
  
  Sigma0 <- rbf_kernel(s0.B, s0.B, par0[2], par0[1]) + par0[3] * diag(length(s0.B))
  Sigma1 <- rbf_kernel(s1.B, s1.B, par1[2], par1[1]) + par1[3] * diag(length(s1.B))
  
  Y0 <- mvrnorm(n.reps, mu = mu0_hat, Sigma = Sigma0)
  Y1 <- mvrnorm(n.reps, mu = mu1_hat, Sigma = Sigma1)
  
  Delta_Bs <- rowMeans(Y1) - rowMeans(Y0)
  p <- mean(Delta_Bs < 0)
  
  list(
    p           = p,
    theta_hat0  = par0[2],
    theta_hat1  = par1[2],
    sigma2_hat0 = par0[1],
    sigma2_hat1 = par1[1],
    v2_hat0     = par0[3],
    v2_hat1     = par1[3],
    beta0_hat   = beta0,
    beta1_hat   = beta1,
    par0        = par0,
    par1        = par1
  )
}


# Fully nonparametric bootstrap
bootstrap_run_procedure <- function(data,
                                    s0.B,
                                    s1.B,
                                    degree = 3,
                                    n.reps = 200,
                                    spline = FALSE,
                                    n.boot = 200) {
  
  # original fit on the full data
  orig_fit <- run_procedure(
    data   = data,
    s0.B   = s0.B,
    s1.B   = s1.B,
    degree = degree,
    n.reps = n.reps,
    spline = spline
  )
  
  # cluster IDs
  unique_studies <- unique(data$study)
  n_studies <- length(unique_studies)
  
  # storage for bootstrap p-hats
  p_boot <- numeric(n.boot)
  
  # optional: progress bar
  pb <- utils::txtProgressBar(min = 0, max = n.boot, style = 3)
  
  for (b in seq_len(n.boot)) {
    # sample studies with replacement
    sampled_ids <- sample(unique_studies, size = n_studies, replace = TRUE)
    
    # rebuild bootstrap dataset with *renumbered* study labels
    boot_pieces <- vector("list", length(sampled_ids))
    for (j in seq_along(sampled_ids)) {
      tmp <- data[data$study == sampled_ids[j], , drop = FALSE]
      tmp$study <- j  # relabel so duplicates stay as separate clusters
      boot_pieces[[j]] <- tmp
    }
    boot_data <- do.call(rbind, boot_pieces)
    
    # run your original procedure on the bootstrap sample
    fit_b <- run_procedure(
      data   = boot_data,
      s0.B   = s0.B,
      s1.B   = s1.B,
      degree = degree,
      n.reps = n.reps,
      spline = spline
    )
    
    p_boot[b] <- fit_b$p
    
    utils::setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # bootstrap variance and SE of p_hat
  var_p <- var(p_boot)
  se_p  <- sqrt(var_p)
  lower_ci <- quantile(p_boot,0.025)
  upper_ci <- quantile(p_boot,0.975)
  
  # return everything
  list(
    orig_fit = orig_fit,  # full original output from run_procedure
    p_boot   = p_boot,    # bootstrap p-hats
    var_p    = var_p,
    se_p     = se_p,
    lower_ci = lower_ci,
    upper_ci = upper_ci
  )
}

# Nonparametric bootstrap in parallel.
bootstrap_run_procedure_parallel <- function(data,
                                             s0.B,
                                             s1.B,
                                             degree = 3,
                                             n.reps = 200,
                                             use_spline = FALSE,
                                             n.boot = 200,
                                             n.cores = detectCores() - 1,
                                             seed = NULL,
                                             knots0, knots1,
                                             boundary_knots0, boundary_knots1) {
  
  unique_studies <- unique(data$study)
  n_studies <- length(unique_studies)
  
  # function that runs ONE bootstrap replicate
  one_boot <- function(b) {
    # (optional) make per-rep seeds reproducible
    # set.seed(seed + b)
    
    sampled_ids <- sample(unique_studies,
                          size = n_studies,
                          replace = TRUE)
    
    # rebuild bootstrap dataset with re-labelled study indices
    boot_pieces <- vector("list", length(sampled_ids))
    for (j in seq_along(sampled_ids)) {
      tmp <- data[data$study == sampled_ids[j], , drop = FALSE]
      tmp$study <- j
      boot_pieces[[j]] <- tmp
    }
    boot_data <- do.call(rbind, boot_pieces)
    
    s0.B.boot <- sample(s0.B, size = length(s0.B), replace = TRUE)
    s1.B.boot <- sample(s1.B, size = length(s1.B), replace = TRUE)
    
    fit_b <- calculate_p_hat(
      data   = boot_data,
      s0.B   = s0.B.boot,
      s1.B   = s1.B.boot,
      degree = degree,
      n.reps = n.reps,
      use_spline = use_spline,
      knots0 = knots0,
      knots1 = knots1, 
      boundary_knots0 = boundary_knots0,
      boundary_knots1 = boundary_knots1
    )
    
    fit_b$p
  }
  
  # run all bootstrap reps in parallel
  p_boot <- unlist(
    mclapply(
      X        = seq_len(n.boot),
      FUN      = one_boot,
      mc.cores = n.cores
    )
  )
  
  var_p <- var(p_boot)
  se_p  <- sqrt(var_p)
  lower_ci <- quantile(p_boot,0.025)
  upper_ci <- quantile(p_boot,0.975)
  
  list(
    p_boot   = p_boot,
    var_p    = var_p,
    se_p     = se_p,
    lower_ci = lower_ci,
    upper_ci = upper_ci
  )
}

# This version has delta + MC over S_B
run_procedure_delta_fixed <- function(
  data,
  s0.B,
  s1.B,
  degree = 3,
  use_spline = FALSE,
  n_mc = 200,
  knots0, knots1, boundary_knots0, boundary_knots1
) {
  
  require(numDeriv)
  require(Matrix)
  if (use_spline) require(splines)
  
  ## ------------------------------------------------------------------
  ## 1. Build study list
  ## ------------------------------------------------------------------
  unique.studies <- unique(data$study)
  studies <- lapply(unique.studies, function(st) {
    list(
      s0.A = data$S[data$study == st & data$G == 0],
      y0.A = data$Y[data$study == st & data$G == 0],
      s1.A = data$S[data$study == st & data$G == 1],
      y1.A = data$Y[data$study == st & data$G == 1]
    )
  })
  
  ## ------------------------------------------------------------------
  ## 3. Negative log-likelihood
  ## ------------------------------------------------------------------
  negloglike <- function(par, studies, treat) {
    sigma2 <- par[1]; theta <- par[2]; v2 <- par[3]
    beta <- par[-c(1:3)]
    
    ll <- 0
    for (st in studies) {
      if (treat) {
        x <- st$s1.A; y <- st$y1.A
      } else {
        x <- st$s0.A; y <- st$y0.A
      }

      X <- make_basis(x,
                      degree = degree,
                      use_spline = use_spline,
                      knots = if(treat) knots1 else knots0,
                      boundary_knots = if(treat) boundary_knots1 else boundary_knots0)
      
      m <- X %*% beta
      C <- rbf_kernel(x, x, theta, sigma2) + v2 * diag(length(x))
      ll <- ll + loglik_contrib(C, y, m)
    }
    ll
  }
  
  ## ------------------------------------------------------------------
  ## 4. Fit parameters
  ## ------------------------------------------------------------------
  p_dim <- if(use_spline) length(knots0) + degree + 1 else degree + 1
  init  <- rep(1, 3 + p_dim)
  
  opt0 <- nlminb(init, negloglike, studies = studies, treat = FALSE)
  opt1 <- nlminb(init, negloglike, studies = studies, treat = TRUE)
  
  par0 <- opt0$par
  par1 <- opt1$par
  
  beta0 <- par0[-c(1:3)]
  beta1 <- par1[-c(1:3)]
  
  ## ------------------------------------------------------------------
  ## 5. B-study mean and covariance
  ## ------------------------------------------------------------------
  X0B <- make_basis(s0.B,
                    degree = degree,
                    use_spline = use_spline,
                    knots = knots0,
                    boundary_knots = boundary_knots0)
  
  X1B <- make_basis(s1.B,
                    degree = degree,
                    use_spline = use_spline,
                    knots = knots1,
                    boundary_knots = boundary_knots1)
  
  
  mu0 <- as.vector(X0B %*% beta0)
  mu1 <- as.vector(X1B %*% beta1)
  
  S0 <- rbf_kernel(s0.B, s0.B, par0[2], par0[1]) + par0[3] * diag(length(s0.B))
  S1 <- rbf_kernel(s1.B, s1.B, par1[2], par1[1]) + par1[3] * diag(length(s1.B))
  
  a0 <- rep(1 / length(s0.B), length(s0.B))
  a1 <- rep(1 / length(s1.B), length(s1.B))
  
  mDelta <- as.numeric(a1 %*% mu1 - a0 %*% mu0)
  vDelta <- as.numeric(t(a1) %*% S1 %*% a1 + t(a0) %*% S0 %*% a0)
  
  p_hat <- pnorm(0, mean = mDelta, sd = sqrt(vDelta))
  
  ## ------------------------------------------------------------------
  ## 6. Delta method
  ## ------------------------------------------------------------------
  H0 <- numDeriv::hessian(function(p) negloglike(p, studies, FALSE), par0)
  H1 <- numDeriv::hessian(function(p) negloglike(p, studies, TRUE), par1)
  
  Sigma_theta <- Matrix::bdiag(solve(H0), solve(H1))
  theta_hat <- c(par0, par1)
  d0 <- length(par0)
  
  mgrad <- function(th) {
    b0 <- th[1:d0][-c(1:3)]
    b1 <- th[(d0 + 1):length(th)][-c(1:3)]
    as.numeric(a1 %*% (X1B %*% b1) - a0 %*% (X0B %*% b0))
  }
  
  grad_m <- numDeriv::grad(mgrad, theta_hat)
  var_mDelta <- as.numeric(t(grad_m) %*% Sigma_theta %*% grad_m)
  
  z <- mDelta / sqrt(vDelta)
  var_p_delta <- (dnorm(z)^2 / vDelta) * var_mDelta
  
  ## ------------------------------------------------------------------
  ## 7. Monte Carlo over S_B
  ## ------------------------------------------------------------------
  p_star <- replicate(n_mc, {
    s0b <- sample(s0.B, replace = TRUE)
    s1b <- sample(s1.B, replace = TRUE)
    
    X0b <- make_basis(s0b,
                      degree = degree,
                      use_spline = use_spline,
                      knots = knots0,
                      boundary_knots = boundary_knots0)
    
    X1b <- make_basis(s1b,
                      degree = degree,
                      use_spline = use_spline,
                      knots = knots1,
                      boundary_knots = boundary_knots1)
    
    mu0b <- X0b %*% beta0
    mu1b <- X1b %*% beta1
    
    S0b <- rbf_kernel(s0b, s0b, par0[2], par0[1]) + par0[3] * diag(length(s0b))
    S1b <- rbf_kernel(s1b, s1b, par1[2], par1[1]) + par1[3] * diag(length(s1b))
    
    vDb <- as.numeric(t(a1) %*% S1b %*% a1 + t(a0) %*% S0b %*% a0)
    pnorm(0, mean(mu1b) - mean(mu0b), sqrt(vDb))
  })
  
  var_sb <- var(p_star)
  
  ## ------------------------------------------------------------------
  ## 8. Return
  ## ------------------------------------------------------------------
  list(
    p_hat = p_hat,
    se_p = sqrt(var_p_delta + var_sb),
    var_delta = var_p_delta,
    var_sb = var_sb,
    mDelta = mDelta,
    vDelta = vDelta,
    par0 = par0,
    par1 = par1,
    lower_ci = p_hat - 1.96*sqrt(var_p_delta + var_sb),
    upper_ci = p_hat + 1.96*sqrt(var_p_delta + var_sb)
    
  )
}


full_procedure <- function(
  data,
  s0.B,
  s1.B,
  degree = 3,
  use_spline = FALSE,
  n.reps = 200,
  calculate_se = TRUE,
  try_analytic = TRUE,
  n_mc = 200,  # This is for the delta + MC
  n_bootstrap = 200 # Fully nonparametric bootstrap
) {
  
  # -------------------------------------------------------------------
  # Set knots if using splines
  # -------------------------------------------------------------------
  if (use_spline) {
    control.vec <- c(s0.B, data$S[data$G == 0])
    treat.vec   <- c(s1.B, data$S[data$G == 1])
    boundary_knots0 <- c(min(control.vec) - 1, max(control.vec) + 1)
    boundary_knots1 <- c(min(treat.vec)   - 1, max(treat.vec)   + 1)
    knots0 <- as.numeric(quantile(control.vec, c(0.33, 0.67)))
    knots1 <- as.numeric(quantile(treat.vec, c(0.33, 0.67)))
  } else {
    knots0 <- knots1 <- boundary_knots0 <- boundary_knots1 <- NULL
  }
  
  # The first thing we do is calculate p_hat.
  p_hat <- calculate_p_hat(data, s0.B, s1.B,
                           degree = degree,
                           n.reps = n.reps,
                           use_spline = use_spline,
                           knots0 = knots0, 
                           knots1 = knots1, 
                           boundary_knots0 = boundary_knots0,
                           boundary_knots1 = boundary_knots1)
  
  if (!calculate_se) {
    return(p_hat$p)
  }
  
  # Now, we try to get the variance using the delta + MC.
  # If it doesn't work, fall back on using nonparametric.
  if (try_analytic) {
    errored_out <- FALSE
    out <- tryCatch(
      {
        tmp <- run_procedure_delta_fixed(data, s0.B, s1.B, degree, use_spline, n_mc,
                                         knots0 = knots0, 
                                         knots1 = knots1, 
                                         boundary_knots0 = boundary_knots0,
                                         boundary_knots1 = boundary_knots1)
        
        if (!is.finite(tmp$se_p) || tmp$se_p <= 0)
          stop("Invalid SE from delta+MC")
        
        tmp
      },
      error = function(e) {
        message("Fallback triggered: ", e$message)
        
        # the parallel bootstrap returns a list with se_p
        errored_out <- TRUE
        bootstrap_run_procedure_parallel(data,
                                         s0.B,
                                         s1.B,
                                         degree,
                                         n.reps,
                                         use_spline,
                                         n_bootstrap,
                                         n.cores = detectCores() - 1,
                                         seed = NULL,
                                         knots0 = knots0, 
                                         knots1 = knots1, 
                                         boundary_knots0 = boundary_knots0,
                                         boundary_knots1 = boundary_knots1)
      }
    )
    return(list(
      p_hat = p_hat$p,
      se_p = out$se_p,
      lower_ci = out$lower_ci,
      upper_ci = out$upper_ci,
      var_method = if(errored_out) "bootstrap" else "analytic"
    ))
  } else {
    out <- bootstrap_run_procedure_parallel(data,
                                            s0.B,
                                            s1.B,
                                            degree,
                                            n.reps,
                                            use_spline,
                                            n_bootstrap,
                                            n.cores = detectCores() - 1,
                                            seed = NULL,
                                            knots0 = knots0, 
                                            knots1 = knots1, 
                                            boundary_knots0 = boundary_knots0,
                                            boundary_knots1 = boundary_knots1)
    return(list(
      p_hat = p_hat$p,
      se_p = out$se_p,
      lower_ci = as.numeric(out$lower_ci),
      upper_ci = as.numeric(out$upper_ci),
      var_method = "bootstrap"
    ))
  }
}


#################################################
## Function that implements Elliott 2015 MLE approach
#################################################

#same input setup as our main function

# data should be a dataframe containing data for all completed studies, with variables: S (the surrogate), Y (the outcome), G (the treatment indicator with 1 as treatment and 0 as control), study (the study indicator)


elliott_prob = function(data, s0.B, s1.B) {
  data$study = as.factor(data$study)
  
  ########################################################
  # Now fitting the Normal model
  ########################################################
  
  # For each trial, estimate treatment effect on S and Y and get sampling covariances
  
  trial_ids = levels(data$study)
  ntrial = length(trial_ids)
  
  # will hold the treatment effects; columns: delta_s, delta_t
  
  theta_hat = matrix(NA, nrow=ntrial, ncol=2)    
  rownames(theta_hat) = trial_ids
  colnames(theta_hat) = c("delta_s","delta_t")
  
  # within-trial sampling covariance matrix for each trial
  
  Sigma_list = vector("list", ntrial)            
  
  for(i in seq_along(trial_ids)){
    tr = trial_ids[i]
    dsub = subset(data, study == tr)
    
    # simple OLS of outcome on G per trial (intercept + G)
    fitS = lm(S ~ G, data = dsub)
    fitY = lm(Y ~ G, data = dsub)
    
    # coefficient on G is the treatment effect estimate
    delta_s_hat = coef(fitS)["G"]
    delta_t_hat = coef(fitY)["G"]
    theta_hat[i, ] = c(delta_s_hat, delta_t_hat)
    
    # Compute empirical pooled variances in treated/control groups 		# and then compute cov across outcomes
    
    z1 = dsub$G == 1
    n1 = sum(z1); n0 = sum(!z1)
    # sample means
    sS1 = var(dsub$S[z1]) 
    sS0 = var(dsub$S[!z1]) 
    sY1 = var(dsub$Y[z1]) 
    sY0 = var(dsub$Y[!z1]) 
    
    # within-group covariance between S and Y
    cov1 = cov(dsub$S[z1], dsub$Y[z1]) 
    cov0 = cov(dsub$S[!z1], dsub$Y[!z1])
    
    var_delta_s = (sS1 / n1) + (sS0 / n0)
    var_delta_t = (sY1 / n1) + (sY0 / n0)
    cov_delta_st = (cov1 / n1) + (cov0 / n0)
    
    Sigma_list[[i]] = matrix(c(var_delta_s, cov_delta_st, 					cov_delta_st, var_delta_t), nrow=2)
  }		
  
  # prepare inputs for mvmeta
  # mvmeta expects a matrix of observed effects and either a 			# list of within-variance matrices or a block-diagonal 				# matrix
  
  Y_obs = theta_hat
  S_within = Sigma_list
  
  # Fit multivariate meta-analytic model via ML (method="ml"). 
  # mvmeta will estimate mu and Tau (between-trial covariance).
  
  fit_mv = mvmeta(Y_obs, S_within, method = "ml")  # ML estimation
  
  # extract estimates
  # Bs, Bt (vector), this matches parameters well
  
  mu_hat = coef(fit_mv) 
  
  # between-trial covariance matrix (2 x 2) 
  
  Tau_hat = fit_mv$Psi            
  
  ########################################################
  # Calculating paradox probability from these estimates
  ########################################################
  
  
  observed.delta.s = mean(s1.B) - mean(s0.B) 
  
  #use properties of Y|X=x for bivariate normal X and Y
  
  rho = Tau_hat[1,2]/(sqrt(Tau_hat[1,1])*sqrt(Tau_hat[2,2]))
  cond.mu =  mu_hat[2] + rho*sqrt(Tau_hat[2,2])*((observed.delta.s-	mu_hat[1])/sqrt(Tau_hat[1,1]))
  cond.var = (1-rho^2)*Tau_hat[2,2]
  
  conditional.prob.paradox = pnorm(0, cond.mu, sqrt(cond.var))
  
  return(conditional.prob.paradox)
  
  
}

