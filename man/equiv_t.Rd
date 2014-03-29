\name{equiv_t}
\alias{equiv_t}
\title{Schuirmann's Test of the Equivalence of Two Independent Groups}
\usage{
equiv_t(x, y, equivint, tr = 0.2, alpha = 0.05, varequal = FALSE, ...)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{a numeric vector}

  \item{equivint}{equivalence interval}

  \item{tr}{proportion of data to trim. When \code{tr == 0}
  the standard Schuirmann test is performed}

  \item{alpha}{desired alpha level}

  \item{varequal}{logical; assume equal variances? Only
  applicable when tr == 0}

  \item{...}{additional arguments to be passed}
}
\value{
returns a \code{data.frame} containing the directional
t-statitics and their associated degrees of freedom and
p-values
}
\description{
Perform Schuirmann's two-sample equivalence test,
potentially with Yuen's formula for trimmed means on the
data in x and y if trimming is included. Missing values are
automatically removed.
}
\examples{
\dontrun{
#equivalence test between v1 and v2 with an interval of .2
v1 <- rnorm(100)
v2 <- v1 + rnorm(100, 2)
equiv_t(v1, v2, .2)
equiv_t(v1, v2, .2, tr=0, varequal=TRUE)
equiv_t(v1, v2, .2, tr=0, varequal=FALSE)
}
}
\author{
Rob Cribbie \email{cribbie@yorku.ca} and Phil Chalmers
\email{rphilip.chalmers@gmail.com}
}

