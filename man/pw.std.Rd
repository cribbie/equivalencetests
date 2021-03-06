% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pw.std.R
\name{pw.std}
\alias{print.pw.std}
\alias{pw.std}
\title{Pairwise standardized}
\usage{
pw.std(data, repeated, ei, alpha = 0.05)

\method{print}{pw.std}(x, ...)
}
\arguments{
\item{data}{a data.frame object}

\item{repeated}{a character vector containing the names of the repeated measures variables}

\item{ei}{equivalence interval, in standardized metric}

\item{alpha}{alpha/significance level}

\item{x}{object of class \code{pw.std}}

\item{...}{additional arguments}
}
\description{
An equivalence test for repeated measures. Tests the null hypothesis that at least one pairwise mean difference 
is not practically equivalent, as defined by an equivalence interval. This null hypothesis is constructed as an intersection-union test, so the null hypothesis fails to be rejected if even one pairwise mean difference is not found to be statistically equivalent.
}
\examples{
dat <- simRepDat()
pw.std(data=dat, repeated=c('repA', 'repB', 'repC'), ei=1)
}

