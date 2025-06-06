\name{monotonicity_test}
\alias{monotonicity_test}

\title{
Monotonicity test
}
\description{
Runs the test of monotonicity for a regression function.
}
\usage{
monotonicity_test(X, Y, h = NA, m = 5, bootstrap_n = 100, 
alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Vector of X values.
}
  \item{Y}{
Vector of Y values.
}
  \item{h}{
Bandwidth for the kernel smoother.
}
  \item{m}{
Window size to calculate linear regression.
}
  \item{bootstrap_n}{
Desired number of bootstrap samples.
}
  \item{alpha}{
Desired alpha level of the test.
}

}
\value{
\item{T_m_value}{Value of the test statistic.}
\item{p_val}{p-value of test}
\item{reject}{whether the test rejects the null}
\item{T_m_samples}{Vector of bootstrapped values of test statistic}
}
\references{
Hall, Peter, and Nancy E. Heckman. "Testing for monotonicity of a regression mean by calibrating for linear functions." Annals of Statistics (2000): 20-39.
}
\author{
Emily Hsiao
}

\keyword{internal}
