\name{nnr_test}
\alias{nnr_test}

\title{
Non-negative residual treatment effect function
}
\description{
Runs the test of non-negative residual treatment effect given two sets of surrogate and primary endpoint values.
}
\usage{
nnr_test(s0, y0, s1, y1, n_bootstrap = 200, alpha = 0.05)
}

\arguments{
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
Vector of endpoint values in treatment group.
}
  \item{n_bootstrap}{
Desired number of bootstrap samples.
}
  \item{alpha}{
Desired alpha level of test.
}
}

\value{
\item{p_value}{p-value of test.}
\item{reject}{Whether the test rejects the null hypothesis.}
\item{s_hat}{Calculated value of test statistic.}
\item{s_vec}{Vector of bootstrapped values of test statistic.}
}
\references{
Hsiao et all 2024 (Under review)
}
\author{
Emily Hsiao
}

\keyword{internal}

