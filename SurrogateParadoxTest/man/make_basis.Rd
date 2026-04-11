\name{make_basis}
\alias{make_basis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Basis
}
\description{
Make basis matrix; used by main meta-analytic function, not intended to be called directly by the user.
}
\usage{
make_basis(x, degree = 3, use_spline = FALSE, knots = NULL, boundary_knots = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of values.
}
  \item{degree}{
Polynomial Degree.
}
  \item{use_spline}{
Use B-spline or monomial.
}
  \item{knots}{
Knots.
}
  \item{boundary_knots}{
Boundary knots.
}
}

\value{
Basis matrix.
}
