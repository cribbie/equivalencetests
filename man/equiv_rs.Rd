\name{equiv_rs}
\alias{equiv_rs}
\title{Test whether correlations are independent}
\usage{
  equiv_rs(dat1, dat2, equiv_int, n1 = NULL, n2 = NULL,
    alpha = 0.05)
}
\arguments{
  \item{dat1}{a matrix or data.frame containing raw data
  used to compute the first correlation.  A scalar input
  may be used as well specifying the sample correlation
  directly}

  \item{dat2}{a matrix or data.frame containing raw data
  used to compute the second correlation.  A scalar input
  may be used as well specifying the sample correlation
  directly}

  \item{equivint}{equivalence interval for correlation}

  \item{alpha}{desired alpha level}

  \item{na.rm}{logical; remove missing values?}

  \item{...}{additional arguments to be passed}
}
\description{
  Given raw data or known sample correlations, test whether
  correlation values are equivalent within a specified
  interval. The null hypothesis is that the groups are not
  equivalent.
}
\details{
  Function uses Anderson and Hauck's (1983) equivalence
  test. Because the p value derived from the test is only
  an approximation, the CIs may produce results that fall
  outside of the equivalence interval at small sample
  sizes, even with a statistically significant p value.
  They therefore provide a measure of precision, but cannot
  be used as a reject/fail to reject decision regarding the
  null hypothesis.
}
\examples{
\dontrun{
#raw data
set.seed(1234)
dat1 <- cbind(rnorm(100), rnorm(100))
dat2 <- cbind(rnorm(200), rnorm(200))
equiv_rs(dat1, dat2, .2)

#correlations input
r1 <- cor(dat1)
r2 <- cor(dat2)
equiv_rs(r1, r2, .2, n1 = 100, n2 = 200)
}
}
\author{
  Rob Cribbie \email{cribbie@yorku.ca} and Phil Chalmers
  \email{rphilip.chalmers@gmail.com}
}

