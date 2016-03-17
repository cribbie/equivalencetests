#' Schuirmann's Test of the Equivalence of Two Independent Groups
#' 
#' Perform Schuirmann's two-sample equivalence test, potentially with Yuen's formula for
#' trimmed means on the data in x and y if trimming is included. Missing values are
#' automatically removed.
#' 
#' @aliases equiv_t
#' @param x a numeric vector
#' @param y a numeric vector
#' @param equivint equivalence interval
#' @param tr proportion of data to trim. When \code{tr == 0} the standard Schuirmann test
#'   is performed
#' @param alpha desired alpha level
#' @param varequal logical; assume equal variances? Only applicable when tr == 0
#' @param ... additional arguments to be passed
#' 
#' @return returns a \code{data.frame} containing the directional t-statitics and 
#'   their associated degrees of freedom and p-values
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export equiv_t
#' @examples
#' \dontrun{
#' #equivalence test between v1 and v2 with an interval of .2
#' v1 <- rnorm(100)
#' v2 <- v1 + rnorm(100, 2)
#' equiv_t(v1, v2, .2)
#' equiv_t(v1, v2, .2, tr=0, varequal=TRUE)
#' equiv_t(v1, v2, .2, tr=0, varequal=FALSE)
#' }
equiv_t <- function(x, y, equivint, tr = 0.2, alpha = 0.05, varequal = FALSE, 
    ...) {
    # Perform Schuirmann's two-sample equivalence test using Yuen's
    # formula for trimmed means on the data in x and y.  The
    # default amount of trimming is 20% Missing values (values
    # stored as NA) are automatically removed.  The p-values are
    # returned in yuen$pvals
    equivt <- function(x, y, equivint, alpha = 0.05, varequal = FALSE, 
        na.rm = TRUE, ...) {
        if (na.rm) 
            x <- x[!is.na(x)]
        if (na.rm) 
            y <- y[!is.na(y)]
        if (varequal == FALSE) {
            t1 <- (mean(x) - mean(y) - equivint)/sqrt((var(x)/length(x)) + 
                (var(y)/length(y)))
            t2 <- (mean(x) - mean(y) + equivint)/sqrt((var(x)/length(x)) + 
                (var(y)/length(y)))
            dft <- (((var(x)/length(x)) + (var(y)/length(y)))^2)/((var(x)^2/(length(x)^2 * 
                (length(x) - 1))) + (var(y)^2/(length(y)^2 * (length(y) - 
                1))))
            probt1 <- pt(t1, dft, lower.tail = T)
            probt2 <- pt(t2, dft, lower.tail = F)
            ifelse(probt1 < alpha & probt2 < alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
                decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
        }
        if (varequal == TRUE) {
            t1 <- (mean(x) - mean(y) - equivint)/sqrt(((((length(x) - 
                1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2))/(length(x) + 
                length(y) - 2)) * (1/length(x) + 1/length(y)))
            t2 <- (mean(x) - mean(y) + equivint)/sqrt(((((length(x) - 
                1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2))/(length(x) + 
                length(y) - 2)) * (1/length(x) + 1/length(y)))
            dft <- length(x) + length(y) - 2
            probt1 <- pt(t1, dft, lower.tail = T)
            probt2 <- pt(t2, dft, lower.tail = F)
            ifelse(probt1 < alpha & probt2 < alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
                decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
        }
        title <- "Schuirmann's Test of the Equivalence of Two Independent Groups"
        ei <- (c(equivint))
        names(ei) <- c("equivalence interval")
        tstats <- c(t1, t2)
        dfs <- c(dft, dft)
        pvals <- c(probt1, probt2)
        names(tstats) <- c("t1", "t2")
        names(dfs) <- c("dft1", "dft2")
        names(pvals) <- c("p_t1", "p_t2")
        return(data.frame(t_ests = tstats, df = dfs, pvals = pvals))
    }
    
    # begin
    if (tr == 0.5) 
        stop("Using tr=.5 is not allowed; use a method designed for medians")
    if (tr > 0.25) 
        print("Warning: with tr>.25 type I error control might be poor")
    x <- x[!is.na(x)]  # Remove any missing values in x
    y <- y[!is.na(y)]  # Remove any missing values in y
    if (tr == 0) 
        return(equivt(x = x, y = y, equivint = equivint, alpha = alpha, 
            varequal = varequal, na.rm = FALSE, ...))
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
    return(data.frame(t_tests = tests, df = df, pvals = pvals))
} 
