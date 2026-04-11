\name{calculate_p_hat}
\alias{calculate_p_hat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate resilience probability.
}
\description{
Calculates resilience probability; used by main meta-analytic function, not intended to be called directly by the user.
}
\usage{
calculate_p_hat(data, s0.B, s1.B, degree = 3, n.reps = 200, use_spline = FALSE, 
knots0, knots1, boundary_knots0, boundary_knots1)
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
Vector of Y values in Study B in treatment group.
}
  \item{degree}{
Degree of basis expansion.
}
  \item{n.reps}{
Number of Delta_B samples.
}
  \item{use_spline}{
Use B-spline basis.
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
\item{p}{Estimated resilience probability.}
\item{theta_hat0}{Estimated theta in control group.}
\item{theta_hat1}{Estimated theta in treatment group.}
\item{sigma2_hat0}{Estimated sigma^2 in control group.}
\item{sigma2_hat1}{Estimated sigma^2 in treatment group.}
\item{v2_hat0}{Estimated v^2 in control group.}
\item{v2_hat1}{Estimated v^2 in treatment group.}
\item{beta0_hat}{Estimated beta vector in control group.}
\item{beta1_hat}{Estimated beta vector in treatment group.}
\item{par0}{Estimated parameter vector for control group.}
\item{par1}{Estimated parameter vector for treatment group.}
}
