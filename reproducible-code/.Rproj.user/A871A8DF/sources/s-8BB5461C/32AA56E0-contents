\name{modified_S_stat}
\alias{modified_S_stat}

\title{
Helper function for non-negative residual treatment effect test.
}
\description{
Calculates the value of the test statistic in the NNR test.
}
\usage{
modified_S_stat(mu0_hat, mu1_hat, s0, y0, s1, y1, grid_x, boot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0_hat}{
Kernel-smoothed estimate of \eqn{\mu_0} function.
}
  \item{mu1_hat}{
Kernel-smoothed estimate of \eqn{\mu_1} function.
}
  \item{s0}{
Vector of surrogate values in control group.
}
  \item{y0}{
Vector of endpoint values in control group.
}
  \item{s1}{
Vector of surrogate values in treatment group.
}
  \item{y1}{
Vector of endpoing values in the treatment group.
}
  \item{grid_x}{
Values of s over which supremum is calculated.
}
  \item{boot}{
Whether this is a bootstrapped statistic or the test statistic. 
}
}

\value{
\item{s_hat}{Value of the test statistic}
\item{sup}{Value of the supremum over grid_x}
}

\author{
Emily Hsiao
}

\keyword{internal}
