#' Test whether correlations are independent
#' 
#' Given raw data or known sample correlations, test whether correlation values
#' are equivalent within a specified interval. The null hypothesis is that the groups
#' are not equivalent.
#' 
#' Function uses Anderson and Hauck's (1983) equivalence test. Because
#' the p value derived from the test is only an approximation, the CIs may produce
#' results that fall outside of the equivalence interval at small sample sizes,
#' even with a statistically significant p value. They therefore provide a measure
#' of precision, but cannot be used as a reject/fail to reject decision regarding
#' the null hypothesis.
#' 
#' @aliases equiv_rs
#' @param dat1 a matrix or data.frame containing raw data used to compute the first correlation. 
#'   A scalar input may be used as well specifying the sample correlation directly
#' @param dat2 a matrix or data.frame containing raw data used to compute the second correlation. 
#'   A scalar input may be used as well specifying the sample correlation directly
#' @param equiv_int equivalence interval for correlation
#' @param n1 sample size for first correlation (required when input is correlation)
#' @param n2 sample size for second correlation (required when input is correlation)
#' @param alpha desired alpha level
#' @param ... additional arguments to be passed
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export equiv_corr
#' @examples
#' \dontrun{
#' #raw data
#' set.seed(1234)
#' dat1 <- cbind(rnorm(100), rnorm(100))
#' dat2 <- cbind(rnorm(200), rnorm(200))
#' equiv_rs(dat1, dat2, .2)
#' 
#' #correlations input
#' r1 <- cor(dat1)
#' r2 <- cor(dat2)
#' equiv_rs(r1, r2, .2, n1 = 100, n2 = 200)
#' }
equiv_rs <- function(dat1, dat2, equiv_int, n1 = NULL, n2 = NULL, alpha = 0.05) {
    if(length(dat1) > 1L){
        dat1 <- na.omit(dat1)
        r1 <- cor(dat1)[1,2]
        n1 <- nrow(dat1)
    } else {
        if(is.null(n1))
            stop('Correlation inputs require their respective sample sizes')
        r1 <- dat1
        if(r1 >= 1 || r1 <= -1) stop('Invalid correlation input')
    }
    if(length(dat2) > 1L){
        dat2 <- na.omit(dat2)
        r2 <- cor(dat2)[1,2]
        n2 <- nrow(dat2)
    } else {
        if(is.null(n2))
            stop('Correlation inputs require their respective sample sizes')
        r2 <- dat2
        if(r2 >= 1 || r2 <= -1) stop('Invalid correlation input')
    }
    ser <- sqrt(((1 - r1^2)^2/(n1 - 2)) + ((1 - r2^2)^2/(n2 - 2)))
    p.value <- pnorm((abs(r1 - r2) - equiv_int)/ser) - pnorm((-abs(r1 - r2) - equiv_int)/ser)
    upper <- (r1 - r2) + qnorm(alpha) * ser
    lower <- (r1 - r2) - qnorm(alpha) * ser
    if (lower < upper) {
        lower2 <- lower
        upper2 <- upper
    }
    if (lower > upper) {
        lower2 <- upper
        upper2 <- lower
    }
    reject <- p.value <= alpha
    ret <- data.frame(r1=r1, r2=r2, equiv_interval=equiv_int, 
                      lowerCI=lower2, upperCI=upper2, p=p.value,
                      reject_equivalence=!reject)
    ret
}
