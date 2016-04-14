#' Test equivalence of two groups'correlations/covariances
#' 
#' Given raw data or known sample correlations, test whether correlation/covariance values
#' are equivalent within a specified interval. The null hypothesis is that the groups
#' are not equivalent in degree of association between two variables.
#' 
#' Function uses Anderson and Hauck's (1983) equivalence test. Because
#' the p-value derived from the test is only an approximation, the CIs may produce
#' results that fall outside of the equivalence interval at small sample sizes,
#' even with a statistically significant p-value. They therefore provide a measure
#' of precision, but cannot be used as a reject/fail to reject decision regarding
#' the null hypothesis.
#' 
#' @aliases equiv_rs
#' @param dat1 a matrix or data.frame containing raw data used to compute the first correlation/
#'   covariance. A scalar input may be used as well specifying the sample correlation directly
#' @param dat2 a matrix or data.frame containing raw data used to compute the second correlation/
#'  covariance. A scalar input may be used as well specifying the sample correlation directly
#' @param equiv_int equivalence interval
#' @param n1 sample size for first covariance set (required when input is a correlation)
#' @param n2 sample size for second covariance set (required when input is a correlation)
#' @param betas logical; compare raw beta regression coefficients rather than correlations?
#' @param alpha desired alpha level
#' @param ... additional arguments to be passed
#' 
#' @return returns a list containing the coefficients, the equivalence interval, p-value, and
#'   statistical decision
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Alyssa Counsell \email{counsela@@yorku.ca}
#' @export equiv_rs
#' @examples
#' \dontrun{
#' #raw data
#' set.seed(1234)
#' dat1 <- cbind(rnorm(100), rnorm(100))
#' dat2 <- cbind(rnorm(200), rnorm(200))
#' equiv_rs(dat1, dat2, .2)
#' equiv_rs(dat1, dat2, .2, betas = TRUE)
#' 
#' par(mfrow=c(1,2))
#' plot(dat1)
#' plot(dat2)
#' par(mfrow=c(1,1))
#' 
#' #correlations input
#' r1 <- cor(dat1)[2,1]
#' r2 <- cor(dat2)[2,1]
#' equiv_rs(r1, r2, .2, n1 = 100, n2 = 200)
#' 
#' }
equiv_rs <- function(dat1, dat2, equiv_int, n1 = NULL, 
    n2 = NULL, betas = FALSE, alpha = 0.05) {
    if (betas && (length(dat1) == 1L || length(dat2) == 
        1L)) 
        stop("beta comparisons require raw data inputs")
    if (length(dat1) > 1L) {
        dat1 <- na.omit(dat1)
        if (betas) {
            mod1 <- lm(dat1[, 2] ~ dat1[, 1])
            r1 <- summary.lm(mod1)$coefficients[2, 
                1]
            se1 <- summary.lm(mod1)$coefficients[2, 
                2]
        } else {
            r1 <- cor(dat1)[1, 2]
        }
        n1 <- nrow(dat1)
    } else {
        if (is.null(n1)) 
            stop("Correlation inputs require their respective sample sizes")
        r1 <- dat1
        if (r1 >= 1 || r1 <= -1) 
            stop("Invalid correlation input")
    }
    if (length(dat2) > 1L) {
        dat2 <- na.omit(dat2)
        if (betas) {
            mod2 <- lm(dat2[, 2] ~ dat2[, 1])
            r2 <- summary.lm(mod2)$coefficients[2, 
                1]
            se2 <- summary.lm(mod2)$coefficients[2, 
                2]
        } else {
            r2 <- cor(dat2)[1, 2]
        }
        n2 <- nrow(dat2)
    } else {
        if (is.null(n2)) 
            stop("Correlation inputs require their respective sample sizes")
        r2 <- dat2
        if (r2 >= 1 || r2 <= -1) 
            stop("Invalid correlation input")
    }
    if (betas) {
        ser <- sqrt(se1^2 + se2^2)
    } else {
        ser <- sqrt(((1 - r1^2)^2/(n1 - 2)) + ((1 - 
            r2^2)^2/(n2 - 2)))
    }
    p.value <- pnorm((abs(r1 - r2) - equiv_int)/ser) - 
        pnorm((-abs(r1 - r2) - equiv_int)/ser)
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
    
    ifelse(p.value <= alpha, check_equiv <- "Reject in favour equivalence", 
        check_equiv <- "Do not reject null hypothesis.")
    if (betas) {
        cfs <- data.frame(b1 = r1, b2 = r2)
    } else {
        cfs <- data.frame(r1 = r1, r2 = r2)
    }
    ret <- cbind(cfs, data.frame(equiv_interval = equiv_int, 
        lowerCI = lower2, upperCI = upper2, pValue = p.value, 
        decision = check_equiv))
    class(ret) <- "equiv_rs"
    return(ret)
}
#' @S3method print equiv_rs
#' @rdname equiv_rs
#' @method print equiv_rs
#' @param x object of class \code{equiv_rs}
print.equiv_rs <- function(x, ...) {
    cat("-----------Test equivalence of two groups correlations/covariances ---------\n\n")
    
    cat("Relation 1 = ", x[[1]], "\n")
    cat("Relation 2 = ", x[[2]], "\n")
    cat("Equivalence interval = ", x[[3]], "\n")
    cat("Lower CI = ", x[[4]], "\n")
    cat("Upper CI = ", x[[5]], "\n")
    cat("p-value = ", (round(x[[6]], 4)), "\n")
    cat("Decision = ", x[[7]], "\n")
} 
