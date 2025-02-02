% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblogreg.R
\name{coef.blblogreg}
\alias{coef.blblogreg}
\title{Obtain the coefficients of the logistic model}
\usage{
\method{coef}{blblogreg}(object, ...)
}
\arguments{
\item{object}{a blblogreg object, the fitted model}

\item{...}{other conditions}
}
\value{
return the of object (the fitted model)
}
\description{
Obtain the coefficients of the logistic model
}
\examples{
coef(blblogreg(Direction ~ Lag1 + Lag2 + Volume, data = ISLR::Smarket, m = 3, B = 100, parallel = TRUE))
}
