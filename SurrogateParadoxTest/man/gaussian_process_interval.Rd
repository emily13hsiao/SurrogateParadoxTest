\name{gaussian_process_interval}
\alias{gaussian_process_interval}

\title{
  Gaussian Process Resilience Interval
}
\description{
   Calculates the resilience probability and resilience bound, 
  assuming that the mean function of Study B data is generated according to a 
  Gaussian Process with specified parameters.
}
\usage{
  gaussian_process_interval(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sigma2, theta, 
  n.iter = 500, M = 100, q_quant = 0.1, plot = FALSE, intervals = TRUE, get_var = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s0.A}{
    Vector of surrogate values in the control group of Study A.
  }
  \item{y0.A}{
    Vector of primary outcome values in the control group of Study A.
  }
  \item{s1.A}{
    Vector of surrogate values in the treatment group of Study A.
  }
  \item{y1.A}{
    Vector of primary outcome values in the treatment group of Study A.
}
  \item{s0.B}{
    Vector of surrogate values in the control group of Study B.
  }
  \item{s1.B}{
    Vector of surrogate values in the treatment group of Study B.
  }
  \item{sigma2}{
    Variance parameter of the Gaussian Process.
  }
  \item{theta}{
    Smoothness parameter of the Gaussian Process.
  }
  \item{n.iter}{
    Number of \eqn{\Delta_B} samples to generate.
  }
  \item{M}{
    Number of bootstrap iterations to estimate the SE of \eqn{\hat{p}_0} and \eqn{\hat{q}_\alpha}.
  }
  \item{q_quant}{
    Desired quantile for the resilience bound. Default is 0.10.
  }
  \item{plot}{
    TRUE or FALSE; include plots of randomly generated functions in results. Default is FALSE.
  }
  \item{intervals}{
    TRUE or FALSE; return the inner 95\% of \eqn{\Delta_B} samples.
    Default is TRUE.
  }
  \item{get_var}{
    TRUE or FALSE; return an estimated variance of the \eqn{\hat{p}_0} and \eqn{\hat{q}_\alpha} parameters.
    Default is TRUE.
  }
}

\details{
  More details can be found in Hsiao, Tian, Parast. "Resilience Measures for the Surrogate Paradox." (2025) Under   Review.
}

\value{
  \item{Delta_hats}{Vector of samples of \eqn{\Delta_B}.}
  \item{Delta_estimate}{Mean of samples of \eqn{\Delta_B}.}
  \item{p_hat}{Estimated value \eqn{\hat{p}_0}, the probability of the surrogate paradox.}
  \item{q_hat}{Estimated value \eqn{\hat{q}_\alpha}}, the resilience bound.
  \item{p_se}{Estimated variance of \eqn{\hat{p}_0}; only returned if get_var = TRUE.}
  \item{q_se}{Estimated variance of \eqn{\hat{q}_\alpha}; only returned if get_var = TRUE.}
  \item{control_plot}{Plot including data points of control group of Study A and functions generated for Study B;         only returned if plot = TRUE.}
  \item{treatment_plot}{Plot including data points of treatment group of Study A and functions generated for             Study B; only returned if plot = TRUE.}
  \item{p_interval}{Inner 95\% of \eqn{\hat{p}_0} samples from bootstrap; only returned if intervals = TRUE.}
  \item{q_interval}{Inner 95\% of \eqn{\hat{q}_\alpha} samples from bootstrap; only returned if intervals = TRUE.}
}
\references{
  Hsiao, Tian, Parast. "Resilience Measures for the Surrogate Paradox." (2025) Under Review.
}
\author{
  Emily Hsiao
}

\examples{
  s0.A <- rnorm(20, 3, sqrt(3))
  y0.A <- sapply(s0.A, function(x) 2 * x - 1)
  s1.A <- rnorm(20, 4, sqrt(3))
  y1.A <- sapply(s1.A, function(x) x + 3)
  s0.B <- rnorm(20, 4.75, sqrt(1))
  s1.B <- rnorm(20, 5.25, sqrt(1))
  result <- gaussian_process_interval(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sigma2 = 0.25, 
    theta = 2, n.iter = 20, M = 20)
  result$p_hat
}



