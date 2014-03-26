#' Wellek-Welch equivalence test
#' 
#' This function computes the Wellek-Welch equivalence test for k-groups.
#' 
#' @aliases wellek_welch
#' @param dv a numeric vector for the continous dependent variable
#' @param group grouping variable that can be coerced to a \code{factor}
#' @param eps numeric value defining the size of the equivalance interval
#' @param alpha the appropriate alpha level
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export wellek_welch
#' @examples
#' \dontrun{
#' 
#' dv <- rnorm(1000)
#' group <- rep(paste0('G',1:4), each=length(dv)/4)
#' wellek_welch(dv, group, eps = 0.5)
#' 
#' }
wellek_welch <- function(dv, group, eps, alpha = 0.05) {
    group <- factor(group)
    size <- table(group)
    ng <- length(size)
    owt <- oneway.test(dv ~ group, var.equal = TRUE)
    f <- owt$statistic
    psisq <- f * ((ng - 1)/mean(size))
    crit_psisq <- ((ng - 1)/mean(size)) * qf(p = alpha, df1 = ng - 1, df2 = sum(size) - ng, 
                                             ncp = mean(size) * eps^2)
    pval <- pf(f, df1 = ng - 1, df2 = owt$parameter[2], ncp = mean(size) * eps^2)
    ret <- data.frame(stat = as.numeric(psisq), df1 = ng - 1, df2 = sum(size) - ng, 
                      p.crit = crit_psisq, p.obs = pval, 
                      equivalent = pval < crit_psisq)        
    ret
} 
