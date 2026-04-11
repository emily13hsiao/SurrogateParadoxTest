
# -------------------------------------------------------------------
# Makes the basis matrix
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Makes the basis matrix
# -------------------------------------------------------------------
make_basis <- function(x, degree = 3, use_spline = FALSE, knots = NULL, boundary_knots = NULL) {
  if (!use_spline) {
    # Simple polynomial of given degree
    cbind(1, stats::poly(x, degree = degree, raw = TRUE))
  } else {
    # Spline with specified knots and boundary knots
    if (is.null(knots)) stop("Knots must be provided for splines")
    if (is.null(boundary_knots)) stop("Boundary knots must be provided for splines")
    splines::bs(x, degree = degree, knots = knots,
                Boundary.knots = boundary_knots, intercept = TRUE)
  }
}

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
      mk <- as.vector(Xk %*% beta)
      
      Ck <- rbf_kernel(xk, xk, theta, sigma2) + par[3] * diag(length(xk))
      loglik <- loglik + loglik_contrib(Ck, yk, mk)
    }
    loglik
  }
  
  # -------------------------------------------------------------------
  # Fit parameters
  # -------------------------------------------------------------------
  p_dim <- if(use_spline) length(knots0) + degree + 1 else degree + 1
  initial_params <- rep(1, 3 + p_dim)
  
  opt0 <- stats::nlminb(start = initial_params,
                        objective = negloglike,
                        studies   = studies,
                        treat     = FALSE,
                        lower     = c(0.1, 0.1, 0.1, rep(-Inf, p_dim)),
                        upper     = rep(Inf, 3 + p_dim))
  
  opt1 <- stats::nlminb(start = initial_params,
                        objective = negloglike,
                        studies   = studies,
                        treat     = TRUE,
                        lower     = c(0.1, 0.1, 0.1, rep(-Inf, p_dim)),
                        upper     = rep(Inf, 3 + p_dim))
  
  par0 <- opt0$par
  par1 <- opt1$par
  beta0 <- par0[-c(1:3)]
  beta1 <- par1[-c(1:3)]
  
  # Compute mu_hat
  X0B <- make_basis(s0.B, degree = degree, use_spline = use_spline,
                    knots = knots0, boundary_knots = boundary_knots0)
  X1B <- make_basis(s1.B, degree = degree, use_spline = use_spline,
                    knots = knots1, boundary_knots = boundary_knots1)
  
  mu0_hat <- as.vector(X0B %*% beta0)
  mu1_hat <- as.vector(X1B %*% beta1)
  
  # Covariance matrices
  Sigma0 <- rbf_kernel(s0.B, s0.B, par0[2], par0[1]) + par0[3] * diag(length(s0.B))
  Sigma1 <- rbf_kernel(s1.B, s1.B, par1[2], par1[1]) + par1[3] * diag(length(s1.B))
  
  # Simulate Y
  Y0 <- MASS::mvrnorm(n.reps, mu = mu0_hat, Sigma = Sigma0)
  Y1 <- MASS::mvrnorm(n.reps, mu = mu1_hat, Sigma = Sigma1)
  
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

# -------------------------------------------------------------------
# Parallel bootstrap
# -------------------------------------------------------------------
bootstrap_run_procedure_parallel <- function(data,
                                             s0.B,
                                             s1.B,
                                             degree = 3,
                                             n.reps = 200,
                                             use_spline = FALSE,
                                             n.boot = 200,
                                             n.cores = 2,
                                             seed = NULL,
                                             knots0, knots1,
                                             boundary_knots0, boundary_knots1) {
  
  unique_studies <- unique(data$study)
  n_studies <- length(unique_studies)
  
  one_boot <- function(b) {
    sampled_ids <- sample(unique_studies,
                          size = n_studies,
                          replace = TRUE)
    
    boot_pieces <- lapply(seq_along(sampled_ids), function(j) {
      tmp <- data[data$study == sampled_ids[j], , drop = FALSE]
      tmp$study <- j
      tmp
    })
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
  
  p_boot <- unlist(
    parallel::mclapply(
      X        = seq_len(n.boot),
      FUN      = one_boot,
      mc.cores = n.cores
    )
  )
  
  var_p <- stats::var(p_boot)
  se_p  <- sqrt(var_p)
  
  list(
    p_boot   = p_boot,
    var_p    = var_p,
    se_p     = se_p
  )
}

run_procedure_pab <- function(data, s0.B, s1.B,
                              degree = 3, use_spline = FALSE, n_mc = 200,
                              knots0, knots1, boundary_knots0, boundary_knots1) {
  
  # Ensure packages exist
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' required")
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package 'Matrix' required")
  if (use_spline && !requireNamespace("splines", quietly = TRUE)) stop("Package 'splines' required")
  
  # Build study list
  unique.studies <- unique(data$study)
  studies <- lapply(unique.studies, function(st) {
    list(
      s0.A = data$S[data$study == st & data$G == 0],
      y0.A = data$Y[data$study == st & data$G == 0],
      s1.A = data$S[data$study == st & data$G == 1],
      y1.A = data$Y[data$study == st & data$G == 1]
    )
  })
  
  # Negative log-likelihood
  negloglike <- function(par, studies, treat) {
    sigma2 <- par[1]; theta <- par[2]; v2 <- par[3]; beta <- par[-c(1:3)]
    ll <- 0
    for (st in studies) {
      x <- if(treat) st$s1.A else st$s0.A
      y <- if(treat) st$y1.A else st$y0.A
      X <- make_basis(x, degree = degree, use_spline = use_spline,
                      knots = if(treat) knots1 else knots0,
                      boundary_knots = if(treat) boundary_knots1 else boundary_knots0)
      m <- X %*% beta
      C <- rbf_kernel(x, x, theta, sigma2) + v2 * diag(length(x))
      ll <- ll + loglik_contrib(C, y, m)
    }
    ll
  }
  
  # Fit parameters
  p_dim <- if(use_spline) length(knots0) + degree + 1 else degree + 1
  init  <- rep(1, 3 + p_dim)
  opt0 <- stats::nlminb(init, negloglike, studies = studies, treat = FALSE)
  opt1 <- stats::nlminb(init, negloglike, studies = studies, treat = TRUE)
  par0 <- opt0$par; par1 <- opt1$par
  beta0 <- par0[-c(1:3)]; beta1 <- par1[-c(1:3)]
  
  # B-study mean and covariance
  X0B <- make_basis(s0.B, degree = degree, use_spline = use_spline,
                    knots = knots0, boundary_knots = boundary_knots0)
  X1B <- make_basis(s1.B, degree = degree, use_spline = use_spline,
                    knots = knots1, boundary_knots = boundary_knots1)
  
  mu0 <- as.vector(X0B %*% beta0)
  mu1 <- as.vector(X1B %*% beta1)
  S0 <- rbf_kernel(s0.B, s0.B, par0[2], par0[1]) + par0[3] * diag(length(s0.B))
  S1 <- rbf_kernel(s1.B, s1.B, par1[2], par1[1]) + par1[3] * diag(length(s1.B))
  
  a0 <- rep(1 / length(s0.B), length(s0.B))
  a1 <- rep(1 / length(s1.B), length(s1.B))
  mDelta <- as.numeric(a1 %*% mu1 - a0 %*% mu0)
  vDelta <- as.numeric(t(a1) %*% S1 %*% a1 + t(a0) %*% S0 %*% a0)
  
  p_hat <- stats::pnorm(0, mean = mDelta, sd = sqrt(vDelta))
  
  # Delta method
  H0 <- numDeriv::hessian(function(p) negloglike(p, studies, FALSE), par0)
  H1 <- numDeriv::hessian(function(p) negloglike(p, studies, TRUE), par1)
  Sigma_theta <- Matrix::bdiag(solve(H0), solve(H1))
  theta_hat <- c(par0, par1)
  d0 <- length(par0)
  mgrad <- function(th) {
    b0 <- th[1:d0][-c(1:3)]
    b1 <- th[(d0+1):length(th)][-c(1:3)]
    as.numeric(a1 %*% (X1B %*% b1) - a0 %*% (X0B %*% b0))
  }
  grad_m <- numDeriv::grad(mgrad, theta_hat)
  var_mDelta <- as.numeric(t(grad_m) %*% Sigma_theta %*% grad_m)
  z <- mDelta / sqrt(vDelta)
  var_p_delta <- (stats::dnorm(z)^2 / vDelta) * var_mDelta
  
  # Monte Carlo over S_B
  p_star <- replicate(n_mc, {
    s0b <- sample(s0.B, replace = TRUE)
    s1b <- sample(s1.B, replace = TRUE)
    X0b <- make_basis(s0b, degree = degree, use_spline = use_spline, knots = knots0, boundary_knots = boundary_knots0)
    X1b <- make_basis(s1b, degree = degree, use_spline = use_spline, knots = knots1, boundary_knots = boundary_knots1)
    mu0b <- X0b %*% beta0
    mu1b <- X1b %*% beta1
    S0b <- rbf_kernel(s0b, s0b, par0[2], par0[1]) + par0[3] * diag(length(s0b))
    S1b <- rbf_kernel(s1b, s1b, par1[2], par1[1]) + par1[3] * diag(length(s1b))
    vDb <- as.numeric(t(a1) %*% S1b %*% a1 + t(a0) %*% S0b %*% a0)
    stats::pnorm(0, mean(mu1b) - mean(mu0b), sqrt(vDb))
  })
  var_sb <- stats::var(p_star)
  
  list(p_hat = p_hat,
       se_p = sqrt(var_p_delta + var_sb),
       var_delta = var_p_delta,
       var_sb = var_sb,
       mDelta = mDelta,
       vDelta = vDelta,
       par0 = par0,
       par1 = par1)
}


meta_analytic_resilience <- function(
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
    return(list(p_hat = p_hat$p))
  }
  
  # Now, we try to get the variance using the delta + MC.
  # If it doesn't work, fall back on using nonparametric.
  if (try_analytic) {
    errored_out <- FALSE
    out <- tryCatch(
      {
        tmp <- run_procedure_pab(data, s0.B, s1.B, degree, use_spline, n_mc,
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
                                         n.cores = 2,
                                         seed = NULL,
                                         knots0 = knots0, 
                                         knots1 = knots1, 
                                         boundary_knots0 = boundary_knots0,
                                         boundary_knots1 = boundary_knots1)
      }
    )
    if(errored_out) {
      conf_p = c(as.numeric(quantile(out$p_boot,0.025)), as.numeric(quantile(out$p_boot,0.975)))}
    if(!errored_out) {
      conf_p = c(p_hat$p - 1.96*out$se_p, p_hat$p + 1.96*out$se_p)}
    return(list(
      p_hat = p_hat$p,
      se_p = out$se_p,
      conf_p = conf_p,
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
                                            n.cores = 2,
                                            seed = NULL,
                                            knots0 = knots0, 
                                            knots1 = knots1, 
                                            boundary_knots0 = boundary_knots0,
                                            boundary_knots1 = boundary_knots1)
    return(list(
      p_hat = p_hat$p,
      se_p = out$se_p,
      conf_p = c(as.numeric(quantile(out$p_boot,0.025)), as.numeric(quantile(out$p_boot,0.975))),
      var_method = "bootstrap"
    ))
  }
}



