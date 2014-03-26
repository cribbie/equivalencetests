\name{yuen_equiv}
\alias{yuen_equiv}
\title{Yuen two-sample equivalance test}
\usage{
yuen_equiv(x, y, equivint, tr = 0.2, alpha = 0.05)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{a numeric vector}

  \item{equivint}{equivalence interval}

  \item{tr}{proportion of data to trim}

  \item{alpha}{desired alpha level}

  \item{...}{additional arguments to be passed}
}
\description{
Perform Schuirmann's two-sample equivalence test using
Yuen's formula for trimmed means on the data in x and y.
Missing values (values stored as NA) are automatically
removed.
}
\examples{
\dontrun{
#equivalence correlation test between v1 and v2 with an interval of .2
v1 <- rnorm(100)
v2 <- v1 + rnorm(100, 2)
yuen_equiv(v1, v2, .2)
}
}
\author{
Rob Cribbie \email{cribbie@yorku.ca} and Phil Chalmers
\email{rphilip.chalmers@gmail.com}
}

