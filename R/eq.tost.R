#' Robust and Non-Robust Variants for the two one-sided equivalence tests (TOST): Test equivalence between 2 independent groups
#' 
#' This R function allows for the computation of the original Schuirmann two-sided t-test (TOST) for equivalence, the modified Schuirmann-Welch test (which does not require the variances to be equal), or the Schuirmann-Yuen test (if normality cannot be assumed), 
#' depending on the arguments specified. The original Schuirmann assumes equal variances and normality. The Schuirmann-Welch assumes normality but not equal variances. The Schuirmann-Yuen accounts for unequal variances and nonnormality using trimmed means and Winsorized variances.
#' @aliases tost
#' 
#' x, y, equivint, varequiv = FALSE, normality = FALSE, tr = 0.2, alpha = 0.05,  na.rm = TRUE, ...
#' 
#' @param x a numeric vector for first sample
#' @param y a numeric vector for the second sample
#' @param equivint numeric value defining the size of the equivalance interval
#' @param varequal logical; If true, equal variances are assumed. Only applicable when tr == 0
#' @param normality logical; If true, normality of x and y are assumed. 
#' @param tr proportion of data to trim. When \code{tr == 0}, equivintther the standard Schuirmann test
#'  or Schuirmann-Welch is performed
#' @param alpha the appropriate alpha level
#' @param ... additional arguments to be passed
#' 
#' @return returns a \code{list} 
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#' @export eq.tost
#' @references 
#' 
#' Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. \emph{Journal of pharmacokinetics and biopharmaceutics}, 15(6), 657-680.
# 
#' @examples
#' \dontrun{
#' 
#' x <- rnorm(100, 1)
#' y <- rnorm(100, 1)

#' # Original Schuirman TOST
#' eq.tost(x, y, equivint=0.5, alpha = 0.05, varequiv=FALSE, normality=FALSE)

#' # Schuirmann-Welch
#' eq.tost(x, y, equivint=0.5, alpha = 0.05, varequiv=FALSE, normality=TRUE)

#' # Schuirmann-Yuen TOST (default)
#' eq.tost(x, y, equivint=0.5, alpha = 0.05)

#' }

eq.tost <- function(x, y, equivint, varequiv = FALSE, 
    normality = FALSE, tr = 0.2, alpha = 0.05, 
    na.rm = TRUE, ...) {
    if (na.rm) {
        x <- x[!is.na(x)]
        y <- y[!is.na(y)]
    }
    if (normality) {
        if (varequiv == FALSE) {
            denom <- sqrt((var(x)/length(x)) + 
                (var(y)/length(y)))
            t1 <- (mean(x) - mean(y) - 
                equivint)/denom
            t2 <- (mean(x) - mean(y) + 
                equivint)/denom
            dft <- (((var(x)/length(x)) + 
                (var(y)/length(y)))^2)/((var(x)^2/(length(x)^2 * 
                (length(x) - 1))) + (var(y)^2/(length(y)^2 * 
                (length(y) - 1))))
            probt1 <- pt(t1, dft, lower.tail = T)
            probt2 <- pt(t2, dft, lower.tail = F)
            ifelse(probt1 <= alpha & 
                probt2 <= alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
                decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
            
            title <- "Schuirmann-Welch Test of the Equivalence of Two Independent Groups"
        }
        if (varequiv == TRUE) {
            denom <- sqrt(((((length(x) - 
                1) * sd(x)^2) + ((length(y) - 
                1) * sd(y)^2))/(length(x) + 
                length(y) - 2)) * (1/length(x) + 
                1/length(y)))
            t1 <- (mean(x) - mean(y) - 
                equivint)/denom
            t2 <- (mean(x) - mean(y) + 
                equivint)/denom
            dft <- length(x) + length(y) - 
                2
            probt1 <- pt(t1, dft, lower.tail = T)
            probt2 <- pt(t2, dft, lower.tail = F)
            ifelse(probt1 <= alpha & 
                probt2 <= alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
                decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
            
            title <- "Schuirmann's Test of the Equivalence of Two Independent Groups"
        }
        
    }
    if (normality == FALSE) {
        h1 <- length(x) - 2 * floor(tr * 
            length(x))
        h2 <- length(y) - 2 * floor(tr * 
            length(y))
        print("Winsorized variances are computed.")
        q1 <- (length(x) - 1) * winvar(x, 
            tr)/(h1 * (h1 - 1))
        q2 <- (length(y) - 1) * winvar(y, 
            tr)/(h2 * (h2 - 1))
        dft <- (q1 + q2)^2/((q1^2/(h1 - 
            1)) + (q2^2/(h2 - 1)))
        crit <- qt(1 - alpha/2, dft)
        dif1 <- mean(x, tr) - mean(y, 
            tr) - equivint
        dif2 <- mean(x, tr) - mean(y, 
            tr) + equivint
        t1 <- dif1/sqrt(q1 + q2)
        t2 <- dif2/sqrt(q1 + q2)
        probt1 <- pt(t1, dft)
        probt2 <- 1 - pt(t2, dft)
        ifelse(probt1 <= alpha & probt2 <= 
            alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
            decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
        title <- "Schuirmann-Yuen Test of the Equivalence of Two Independent Groups"
    }
    
    means <- c(mean(x), mean(y))
    names(means) <- c("Mean Grp 1", "Mean Grp 2")
    trimmeans <- c(mean(x, tr), mean(y, 
        tr))
    names(trimmeans) <- c("Trimmed Mean Grp 1", 
        "Trimmed Mean Grp 2")
    sds <- c(sd(x), sd(y))
    names(sds) <- c("SD Grp 1", "SD Grp 2")
    equivint <- (c(equivint))
    names(equivint) <- c("equivalence interval")
    tstats <- c(t1, t2)
    names(tstats) <- c("t1", "t2")
    dfs <- c(dft, dft)
    names(dfs) <- c("dft1", "dft2")
    pvals <- c(probt1, probt2)
    names(pvals) <- c("p_t1", "p_t2")
    res <- list(title, means = means, 
        trimmeans = trimmeans, sds = sds, 
        equivint = equivint, tstats = tstats, 
        dfs = dfs, pvals = pvals, decis = decis)
    return(res)
}
#' @S3method print eq.1way.ww
#' @rdname eq.1way.ww
#' @method print eq.1way.ww
#' @param x object of class \code{eq.1way.ww}
print.eq.1way.ww <- function(x, ...) {
    cat("----", x$title, "----", "\n\n")
    cat("Means:", x$means, "\n")
    cat("SDs:", x$sds, "\n")
    cat("Trimmed Means:", x$trimmeans, 
        "\n")
    cat("The equivalence interval was ", 
        x$equivint, "in unstandardized metric.")
    cat("Test statistics: ", x$tstats, 
        "\n")
    cat("Degrees of freedom: ", x$dfs, 
        "\n")
    cat("p-value = ", x$p.vals, "\n")
    cat("Decision:", x$decis, "\n")
} 
