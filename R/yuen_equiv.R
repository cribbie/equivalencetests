#' Schuirmann/Yuen two-sample equivalance test
#' 
#' Perform Schuirmann's two-sample equivalence test, potentially with Yuen's formula for
#' trimmed means on the data in x and y if trimming is included. Missing values (values stored as NA) are 
#' automatically removed.
#' 
#' @aliases yuen_equiv
#' @param x a numeric vector
#' @param y a numeric vector
#' @param equivint equivalence interval
#' @param tr proportion of data to trim
#' @param alpha desired alpha level
#' @param ... additional arguments to be passed
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export yuen_equiv
#' @examples
#' \dontrun{
#' #equivalence correlation test between v1 and v2 with an interval of .2
#' v1 <- rnorm(100)
#' v2 <- v1 + rnorm(100, 2)
#' yuen_equiv(v1, v2, .2)
#' }
yuen_equiv <- function(x, y, equivint, tr = 0.2, alpha = 0.05) {
    # Perform Schuirmann's two-sample equivalence test using Yuen's formula for trimmed means on the data in
    # x and y.  The default amount of trimming is 20% Missing values (values stored as NA) are automatically
    # removed.  The p-values are returned in yuen$pvals
    if (tr == 0.5) 
        stop("Using tr=.5 is not allowed; use a method designed for medians")
    if (tr > 0.25) 
        print("Warning: with tr>.25 type I error control might be poor")
    x <- x[!is.na(x)]  # Remove any missing values in x
    y <- y[!is.na(y)]  # Remove any missing values in y
    h1 <- length(x) - 2 * floor(tr * length(x))
    h2 <- length(y) - 2 * floor(tr * length(y))
    q1 <- (length(x) - 1) * winvar(x, tr)/(h1 * (h1 - 1))
    q2 <- (length(y) - 1) * winvar(y, tr)/(h2 * (h2 - 1))
    df <- (q1 + q2)^2/((q1^2/(h1 - 1)) + (q2^2/(h2 - 1)))
    crit <- qt(1 - alpha/2, df)
    dif1 <- mean(x, tr) - mean(y, tr) - equivint
    dif2 <- mean(x, tr) - mean(y, tr) + equivint
    test1 <- dif1/sqrt(q1 + q2)
    test2 <- dif2/sqrt(q1 + q2)
    yuenp1 <- pt(test1, df)
    yuenp2 <- 1 - pt(test2, df)
    pvals <- c(yuenp1, yuenp2)
    tests <- c(test1, test2)
    list(tests = tests, pvals = pvals)
} 
