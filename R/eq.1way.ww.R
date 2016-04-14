#' Wellek-Welch equivalence test for k-groups 
#' 
#' This function computes the Wellek-Welch equivalence test for k-groups. This addresses the research question about equivalence for means among k independent groups.
#' The Welch adjustment provides a robust correction for unequal vairances among the k groups.
#' 
#' @aliases eq.1way.ww
#' @param dv a numeric vector for the continous dependent variable
#' @param group grouping variable that can be coerced to a \code{factor}
#' @param eps numeric value defining the size of the equivalance interval
#' @param alpha the appropriate alpha level
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export eq.1way.ww
#' @examples
#' \dontrun{
#' 
#' dv <- rnorm(1000)
#' group <- rep(paste0('G',1:4), each=length(dv)/4)
#' eq.1way.ww(dv, group, eps = 0.5)
#' 
#' }
eq.1way.ww <- function(dv, group, eps, alpha = 0.05) {
    group <- factor(group)
    size <- table(group)
    ng <- length(size)
    owt <- oneway.test(dv ~ group, var.equal = TRUE)
    f <- owt$statistic
    psisq <- f * ((ng - 1)/mean(size))
    crit_psisq <- ((ng - 1)/mean(size)) * qf(p = alpha, 
        df1 = ng - 1, df2 = sum(size) - ng, ncp = mean(size) * 
            eps^2)
    pval <- pf(f, df1 = ng - 1, df2 = owt$parameter[2], 
        ncp = mean(size) * eps^2)
    ifelse(pval < crit_psisq, check_equiv <- "The null hypothesis is rejected.", 
        check_equiv <- "The null hypothesis is not rejected.")
    ret <- data.frame(stat = as.numeric(psisq), df1 = ng - 
        1, df2 = sum(size) - ng, p.crit = crit_psisq, 
        p.obs = pval, decision = check_equiv)
    class(ret) <- "eq.1way.ww"
    return(ret)
}

#' @S3method print eq.1way.ww
#' @rdname eq.1way.ww
#' @method print eq.1way.ww
#' @param x object of class \code{eq.1way.ww}
print.eq.1way.ww <- function(x, ...) {
    cat("----Wellek-Welch equivalence test for k-groups----\n\n")
    cat("F statistic = ", x$stat, "\n")
    cat("with degrees of freedom ", x$df1, ",", x$df2, 
        ".\n")
    cat("Critical value = ", x$p.crit, "\n")
    cat("p-value = ", x$p.obs, "\n")
} 
