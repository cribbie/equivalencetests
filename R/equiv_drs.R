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
#' @param equiv_int equivalence interval
#' @param n sample size when dat input is a vector of correlations
#' @param alpha desired alpha level
#' 
#' @return returns a list containing the p-value, confidence interval, and statistical decision
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export equiv_drs
#' @examples
#' \dontrun{
#' #raw data
#' set.seed(1234)
#' dat <- cbind(rnorm(100), rnorm(100), rnorm(100)) 
#' equiv_drs(dat, equiv_int = .2)
#' 
#' #correlations input
#' r12 <- cor(dat)[2,1]
#' r13 <- cor(dat)[2,1]
#' r23 <- 
#' equiv_drs(c(r12, r13, r23), equiv_int = .2, n = nrow(dat))
#' }
equiv_drs <- function(dat, equiv_int, n = NULL, alpha = 0.05) {
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
    p1 <- pnorm((abs(r12 - r13) - equiv_int) * (sqrt(((n - 1) * (1 + r23))/((2 * ((n - 1)/(n - 3)) * detR) + 
        (((r12 + r13)^2)/4) * ((1 - r23)^3)))))
    p2 <- pnorm((-abs(r12 - r13) - equiv_int) * (sqrt(((n - 1) * (1 + r23))/((2 * ((n - 1)/(n - 3)) * detR) + 
        (((r12 + r13)^2)/4) * ((1 - r23)^3)))))
    p.value <- p1 - p2
    ser <- 1/sqrt(((n - 1) * (1 + r23))/((2 * ((n - 1)/(n - 3)) * detR) + (((r12 + r13)^2)/4) * ((1 - r23)^3)))
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
    out <- list(p.value, CI, decision)
    names(out) <- c("P value", "Confidence Interval of the Difference in Correlations", "Decision")
    print(out)
} 
