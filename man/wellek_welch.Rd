\name{wellek_welch}
\alias{wellek_welch}
\title{Wellek-Welch equivalence test}
\usage{
wellek_welch(dv, group, eps, alpha = 0.05)
}
\arguments{
  \item{dv}{a numeric vector for the continous dependent
  variable}

  \item{group}{grouping variable that can be coerced to a
  \code{factor}}

  \item{eps}{numeric value defining the size of the
  equivalance interval}

  \item{alpha}{the appropriate alpha level}
}
\description{
This function computes the Wellek-Welch equivalence test
for k-groups.
}
\examples{
\dontrun{

dv <- rnorm(100)
group <- rep(paste0('G',1:4), each=25)
wellek_welch(dv, group, eps = 0.25)

}
}
\author{
Rob Cribbie \email{cribbie@yorku.ca} and Phil Chalmers
\email{rphilip.chalmers@gmail.com}
}

