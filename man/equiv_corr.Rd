\name{equiv_corr}
\alias{equiv_corr}
\alias{print.equiv_corr}
\title{Run several equivalence based correlation tests}
\usage{
  equiv_corr(x, y, equivint, alpha = 0.05, na.rm = TRUE,
    ...)

  \method{print}{equiv_corr} (x)
}
\arguments{
  \item{x}{a numeric variable to be correlated with
  \code{y}}

  \item{y}{a numeric variable to be correlated with
  \code{x}}

  \item{equivint}{equivalence interval}

  \item{alpha}{desired alpha level}

  \item{na.rm}{logical; remove missing values?}

  \item{...}{additional arguments to be passed}

  \item{x}{object of class \code{equiv_corr}}
}
\description{
  The following function was created to run a traditional
  Pearson correlation test, an equivalence based test of
  lack of association, an equivalence based test of lack of
  association using a Fisher's z transformation, and an
  equivalence based test of lack of association with
  resampling.
}
\examples{
\dontrun{
#equivalence correlation test between v1 and v2 with an interval of .2
v1 <- rnorm(100)
v2 <- v1 + rnorm(100, 2)
equiv_corr(v1, v2, .2)
}
}
\author{
  Rob Cribbie \email{cribbie@yorku.ca}
}

