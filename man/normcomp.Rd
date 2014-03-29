\name{normcomp}
\alias{normcomp}
\alias{print.normcomp}
\title{Normative comparison tests for evaluating clinical significance}
\usage{
normcomp(premean, presd, pren, postmean, postsd, postn, normmean, normsd, normn,
  alpha = 0.05, ...)

\method{print}{normcomp}(x, ...)
}
\arguments{
  \item{premean}{pretest mean}

  \item{presd}{pretest standard deviation}

  \item{pren}{pretest sample size}

  \item{postmean}{posttest mean}

  \item{postsd}{post test standard deviation}

  \item{postn}{post test sample size}

  \item{normmean}{normal comparison mean}

  \item{normsd}{normal comparison standard deviation}

  \item{normn}{normal comparison sample size}

  \item{alpha}{desired alpha level}

  \item{...}{additional arguments to be passed}

  \item{x}{object of class \code{normcomp}}
}
\description{
The following function conducts normative comparison tests
for evaluating clinical significance, using both the
Cribbie and Arpin-Cribbie (2009) and Kendall et al. (1999)
methods. Normative comparison tests compare a treated
sample to a normal comparison sample to determine if they
are statistically equivalent.
}
\examples{
\dontrun{
#In this example the pretest mean was 12, pretest sd was 4, pretest sample size was 100, posttest mean was 8,
#posttest sd was 5, posttest sample size was 95, normal comparison mean was 7, normal comparison sd was 3,
#normal comparison sample size was 500.
#The desired alpha level is .05, which is the default, so no alpha level is specified
normcomp(12,4,100,8,5,95,7,3,500)
}
}
\author{
Rob Cribbie \email{cribbie@yorku.ca} and Phil Chalmers
\email{rphilip.chalmers@gmail.com}
}
\references{
Cribbie, R. A. & Arpin-Cribbie, C. A. (2009). Evaluating
clinical significance through equivalence testing:
Extending the normative comparisons approach.
\emph{Psychotherapy Research, 19}, 677-686.

Kendall, P. C., Marrs-Garcia, A., Nath, S. R., & Sheldrick,
R. C. (1999). Normative comparisons for the evaluation of
clinical significance. \emph{Journal of Consulting and
Clinical Psychology, 67}, 285-299.
}

