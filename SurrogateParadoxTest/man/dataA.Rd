\name{dataA}
\alias{dataA}
\docType{data}
\title{
Example dataset for meta-analytic analysis
}
\description{
Example dataset for meta-analytic analysis; this dataset contains the surrogate and primary outcome in the treatment and control groups, measured for multiple studies.
}
\usage{data("dataA")}
\format{
  A data frame with 500 observations on the following 4 variables.
  \describe{
    \item{\code{S}}{the surrogate value}
    \item{\code{Y}}{the primary outcome}
    \item{\code{G}}{the treatment indicator}
    \item{\code{study}}{the Study ID}
  }
}

\examples{
data(dataA)
names(dataA)
head(dataA)
}
\keyword{datasets}
