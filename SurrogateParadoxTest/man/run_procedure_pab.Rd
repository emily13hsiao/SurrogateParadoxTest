\name{run_procedure_pab}
\alias{run_procedure_pab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meta-Analytic Resilience Probability using PAB
}
\description{
Calculates the meta-analytic resilience probability, estimating SE using 
partially analytic bootstrap; used by main meta-analytic function, not intended to be called directly by the user.
}
\usage{
run_procedure_pab(data, s0.B, s1.B, degree = 3, use_spline = FALSE, n_mc = 200, 
knots0, knots1, boundary_knots0, boundary_knots1)
}

\arguments{
  \item{data}{
Study A data.
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
  \item{use_spline}{
Whether to use B-spline basis.
}
  \item{n_mc}{
Number of Monte-Carlo iterations.
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
\item{se_p}{Estimated SE of resilience probability.}
\item{var_delta}{Estimated component of SE from Delta method.}
\item{var_sb}{Estimated component of SE from resampling S.}
\item{mDelta}{Mean for Delta method.}
\item{vDelta}{Variance for Delta method.}
\item{par0}{MLE parameters for control group.}
\item{par1}{MLE parameters for treatment group.}

}

