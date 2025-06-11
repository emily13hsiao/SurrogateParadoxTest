\name{gp_resilience_set}
\alias{gp_resilience_set}

\title{
  Gaussian Process Resilience Set
}
\description{
  Creates a plot of the resilience set i.e., the possible parameters of the Gaussian Process such that the 
  probability of the surrogate paradox is below a threshold \eqn{\alpha}.
}
\usage{
  gp_resilience_set(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sigma2_vals, theta_vals, 
  alpha = 0.05)
}

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
  \item{sigma2_vals}{
    Vector of values of \eqn{\sigma^2} for which the probability of surrogate paradox should be calculated with       the given data.
  }
  \item{theta_vals}{
    Vector of values of \eqn{\theta} for which the probability of surrogate paradox should be calculated with         the given data.
  }
  \item{alpha}{
    Threshold value of the surrogate paradox.
  }
}

\details{
  More details can be found in Hsiao, Tian, Parast. "Resilience Measures for the Surrogate Paradox." (2025) Under   Review.
}

\value{
  Returns a ggplot object with values of \eqn{\theta} on the x-axis and \eqn{\sigma^2} on the y-axis, and regions
  where \eqn{P(\Delta_B) < \alpha} with given data and parameters highlighted in blue.
}

\references{
  Hsiao, Tian, Parast. "Resilience Measures for the Surrogate Paradox." (2025) Under Review.
}

\author{
  Emily Hsiao
}

\examples{
n_A=200
n_B=200
s0.A <- rnorm(n_A, mean = 10, sd = 2)
y0.A <- 5 + 0.5 * s0.A + rnorm(n_A, mean = 0, sd = 0.5)

s1.A <- rnorm(n_A, mean = 12, sd = 2.5)
y1.A <- 5.5 + 0.6 * s1.A + rnorm(n_A, mean = 0, sd = 0.6)

 s0.B <- rnorm(n_B, 10, 2)
 s1.B <- rnorm(n_B, 11, 3)

  
  sigma2_values <- seq(0.01, 10, length.out = 30)
  theta_values <- seq(0.01, 10, length.out = 30)
  \donttest{
  gp_resilience_set(s0.A, y0.A, s1.A, y1.A, s0.B, s1.B, sigma2_values, theta_values, 
  alpha = 0.05)
  }
}



