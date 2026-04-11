\name{dataB}
\alias{dataB}
\docType{data}
\title{
Example dataset for meta-analytic analysis, new study
}
\description{
Example dataset for meta-analytic analysis; this dataset contains the surrogate in the treatment (s1.B) and control group (s0.B) for the new study.
}
\usage{data("dataB")}
\format{
 A list of 10 observations in each treatment group:
  \describe{
    \item{\code{s0.B}}{the surrogate values in the control group}
    \item{\code{s1.B}}{the surrogate values in the treated group}
  }}

\examples{
data(dataB)
names(dataB)
head(dataB)
}
\keyword{datasets}
