% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pw.unstd.R
\name{pw.unstd}
\alias{print.pw.unstd}
\alias{pw.unstd}
\title{Pairwise unstandardized}
\usage{
pw.unstd(data, repeated, ei, alpha = 0.05)

\method{print}{pw.unstd}(x, ...)
}
\arguments{
\item{data}{a data.frame object}

\item{repeated}{a character vector containing the names of the repeated measures variables}

\item{ei}{equivalence interval, in unstandardized metric}

\item{alpha}{alpha/significance level}

\item{x}{object of class \code{pw.unstd}}

\item{...}{additional arguments}
}
\description{
Equivalence interval is in unstandardized metric. If at least one pairwise mean difference is found not statistically significant (in equivalence), then the whole test fails to reject the null hypothesis.
}
\examples{
dat <- simRepDat()
pw.unstd(data=dat, repeated=c('repA', 'repB', 'repC'), ei=0.25)
}
\references{
Mara, C. A., & Cribbie, R. A. (2012). Paired-samples tests of equivalence. \emph{Communications in Statistics-Simulation and Computation}, 41(10), 1928-1943.

Wellek, S. (2010). \emph{Testing statistical hypotheses of equivalence and noninferiority}. CRC Press.
}

