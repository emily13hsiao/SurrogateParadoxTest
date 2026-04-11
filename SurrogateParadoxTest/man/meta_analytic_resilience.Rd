\name{meta_analytic_resilience}
\alias{meta_analytic_resilience}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-Analytic Resilience Probability
}
\description{
This function calculates the resilience probability for the meta-analytic setting.
}
\usage{
meta_analytic_resilience(data, s0.B, s1.B, degree = 3, use_spline = FALSE, 
n.reps = 200, calculate_se = TRUE, try_analytic = TRUE, n_mc = 200, n_bootstrap = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data from Study A. Needs to be in the specific format of a table with the columns: 
S, Y, G (treatment indicator), study (study index).
}
  \item{s0.B}{
Vector of surrogate values in Study B in control group.
}
  \item{s1.B}{
Vector of surrogate values in Study B in treatment group.
}
  \item{degree}{
The degree of polynomial desired for the basis expansion.
}
  \item{use_spline}{
TRUE or FALSE: whether to use a B-spline basis expansion for mean function.
}
  \item{n.reps}{
How many samples of Delta_K+1 are desired.
}
  \item{calculate_se}{
TRUE or FALSE: return SE of estimated p. 
}
  \item{try_analytic}{
TRUE or FALSE: run PAB method or default to fully nonparametric bootstrap.
}
  \item{n_mc}{
Number of samples in PAB method.
}
  \item{n_bootstrap}{
Number of samples in nonparametric bootstrap.
}
}
\details{
More details can be found in Hsiao, Parast. 
"A Functional-Class Meta-Analytic Framework for Quantifying Surrogate Resilience" (2026).
}
\value{
  \item{p_hat}{Estimated resilience probability.}
  \item{se_p}{Estimated SE of resilience probability; only if calculate_se = TRUE.}
  \item{conf_p}{95\% confidence interval for resilience probability; only if calculate_se = TRUE.}
  \item{var_method}{Method used to estimate SE, bootstrap or analytic; only if calculate_se = TRUE.}
}

\examples{

#use example datasets
data(dataA)
names(dataA)
head(dataA)

data(dataB)
names(dataB)
head(dataB)

\donttest{
#apply function, no uncertainty quantification
set.seed(1)
result = meta_analytic_resilience(data=dataA, s0.B = dataB$s0.B, s1.B = dataB$s1.B, 
degree = 1, use_spline = FALSE, calculate_se = FALSE)
result

#apply function, uncertainty quantification, analytic SE
set.seed(1)
result = meta_analytic_resilience(data=dataA, s0.B = dataB$s0.B, s1.B = dataB$s1.B, 
degree = 1, use_spline = FALSE, calculate_se = TRUE, try_analytic = TRUE)
result

#apply function, uncertainty quantification, bootstrap SE
set.seed(1)
result = meta_analytic_resilience(data=dataA, s0.B = dataB$s0.B, s1.B = dataB$s1.B, 
degree = 1, use_spline = FALSE, calculate_se = TRUE, try_analytic = FALSE)
result
}
}

