#' Dependent Correlations
#' 
#' Tests the equivalence of the correlations r12 and r13 (variable 1 is overlapping variable) in
#' a three variable correlation matrix. Standard errors based on William's modification 
#' to Hotellings test comparing dependent overlapping correlations.
#' 
#' @aliases equiv_drs
#' @param dat an N x 3 matrix or data.frame containing raw data used to compute the 
#'   correlation matrix between variables. The input may also be a 1 x 3 vector of correlations
#'   (r12, r13, and r23, respectively) and requires a sample size input (N)
#' @param ei equivalence interval
#' @param n sample size when dat input is a vector of correlations
#' @param alpha desired alpha level
#' 
#' @return returns a list containing the p-value, confidence interval, and statistical decision
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Alyssa Counsell \email{counsela@@yorku.ca}
#' @export equiv_drs
#' @examples
#' \dontrun{
#' #raw data
#' set.seed(1234)
#' dat <- cbind(rnorm(100), rnorm(100), rnorm(100)) 
#' equiv_drs(dat, ei = .2)
#' 
#' #correlations input
#' r12 <- cor(dat)[2,1]
#' r13 <- cor(dat)[3,1]
#' r23 <- cor(dat)[3,2]
#' equiv_drs(c(r12, r13, r23), ei = .2, n = nrow(dat))
#' }
equiv_drs <- function(dat, ei, n = NULL, alpha = 0.05) {
    if (length(dat) > 3L) {
        x1 <- dat[, 1L]
        x2 <- dat[, 2L]
        x3 <- dat[, 3L]
        r12 <- cor(x1, x2)
        r13 <- cor(x1, x3)
        r23 <- cor(x2, x3)
        n <- nrow(dat)
    } else {
        r12 <- dat[1L]
        r13 <- dat[2L]
        r23 <- dat[3L]
        if (is.null(n)) 
            stop("sample size input required")
    }
    detR <- (1 - r12^2 - r13^2 - r23^2) + (2 * r12 * r13 * r23)
    p1 <- pnorm((abs(r12 - r13) - ei) * (sqrt(((n - 1) * (1 + r23))/((2 * 
        ((n - 1)/(n - 3)) * detR) + (((r12 + r13)^2)/4) * ((1 - r23)^3)))))
    p2 <- pnorm((-abs(r12 - r13) - ei) * (sqrt(((n - 1) * (1 + r23))/((2 * 
        ((n - 1)/(n - 3)) * detR) + (((r12 + r13)^2)/4) * ((1 - r23)^3)))))
    p.value <- p1 - p2
    ser <- 1/sqrt(((n - 1) * (1 + r23))/((2 * ((n - 1)/(n - 3)) * detR) + 
        (((r12 + r13)^2)/4) * ((1 - r23)^3)))
    upper <- (r12 - r13) + qnorm(alpha) * ser
    lower <- (r12 - r13) - qnorm(alpha) * ser
    if (lower < upper) {
        lower2 <- lower
        upper2 <- upper
    }
    if (lower > upper) {
        lower2 <- upper
        upper2 <- lower
    }
    CI <- c(lower2, upper2)
    decision <- if (p.value <= alpha) {
        "Dependent correlation coefficients can be considered equivalent at alpha"
    } else {
        "Dependent correlation coefficients can NOT be considered equivalent at alpha."
    }
    out <- list(p.value, CI, decision, ei)
    names(out) <- c("P value", "Confidence Interval of the Difference in Correlations", 
        "Decision", "ei")
    class(out) <- "equiv_drs"
    return(out)
}
#' @S3method print equiv_drs
#' @rdname equiv_drs
#' @method print equiv_drs
#' @param x object of class \code{equiv_drs}
print.equiv_drs <- function(x, ...) {
    cat("------Equivalence tests of dependent correlations------\n\n")
    cat("p-value = ", x[[1]], "\n")
    cat("CI = ", x[[2]], "\n")
    cat("Decision: ", x[[3]])
    cat("\n\n")
    cat("Equivalence interval = ", x[[4]])
} 
