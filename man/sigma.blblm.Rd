% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{sigma.blblm}
\alias{sigma.blblm}
\title{Compute the residual standard deviation (sigma) of a fitted model}
\usage{
\method{sigma}{blblm}(object, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{a blblm object, the fitted model}

\item{confidence}{boolean, specify whether to compute the confidence interval}

\item{level}{double, the significance level 1-alpha, 0.95 by default}

\item{...}{other conditions}
}
\value{
return the residual standard deviation (sigma) of the object (the fitted model)
}
\description{
Compute the residual standard deviation (sigma) of a fitted model
}
\examples{
sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
sigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100), confidence = TRUE)
}
