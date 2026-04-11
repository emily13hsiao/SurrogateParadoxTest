\name{bootstrap_run_procedure_parallel}
\alias{bootstrap_run_procedure_parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-Analytic Resilience Probability using Bootstrap
}
\description{
Calculates the meta-analytic resilience probability, estimating the standard error (SE) using the nonparametric bootstrap; used by main meta-analytic function, not intended to be called directly by the user.
}
\usage{
bootstrap_run_procedure_parallel(data, s0.B, s1.B, degree = 3, 
n.reps = 200, 
use_spline = FALSE, 
n.boot = 200, 
n.cores = 2, 
seed = NULL, knots0, knots1, boundary_knots0, boundary_knots1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Dataset in same format as required by meta_analytic_resilience function.
}
  \item{s0.B}{
Vector of S values in Study B in control group.
}
  \item{s1.B}{
Vector of S values in Study B in treatment group.
}
  \item{degree}{
Desired polynomial degree for mean function.
}
  \item{n.reps}{
Number of Delta_B samples to generate.
}
  \item{use_spline}{
Whether to use B-spline basis.
}
  \item{n.boot}{
Number of nonparametric bootstrap samples.
}
  \item{n.cores}{
Number of cores to use in parallelization.
}
  \item{seed}{
Seed for parallelization.
}
  \item{knots0}{
Knots for B-spline in control group.
}
  \item{knots1}{
Knots for B-spline in treatment group.
}
  \item{boundary_knots0}{
Boundary knots for spline in control group.
}
  \item{boundary_knots1}{
Boundary knots for spline in treatment group.
}
}

\value{
\item{p_boot}{Estimated resilience probability.}
\item{var_p}{Estimated variance of resilience probability.}
\item{se_p}{Estimated SE of resilience probability.}
}

