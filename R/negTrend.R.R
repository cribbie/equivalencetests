#' Test of negligible trend (two one-sided tests approach)
#' 
#' This function takes as input the estimates obtained from a linear mixed model. Specifically, some fixed effect slope for time or some other ordered variable, its standard error, and the approximated degrees of freedom are used. 
#' This is essentially a TOST; instead of the mean difference as the parameter of interest, the slope is. 
#' The equivalence interval is then specified with the lower and upper limit bounds for what is considered to be minimum average slope deemed practically important. 
#' For analyses that include random slopes, one would used the fixed effect slope for the test of negligible trend.
#' One could also then evaluate the degree to which individual-level slope coefficients fall within negligible bounds, for assessing individual-level stability. 
#' Eg. What proportion of individual rates of (lack of) change had values that fell within the equivalence interval? 

#' @param  ei equivalence interval
#' @param slope the estimated fixed effect slope obtained from some model
#' @param se the standard error of the slope
#' @param df the degrees of freedom obtained from the fitted model  
#' @export negTrend
negTrend <- function(ei, slope, se, df, n, alpha = 0.05) {
    t1 <- (slope - ei)/se
    t2 <- (slope + ei)/se
    # dft <- summary(mod)$tTable[2,3] #this df is
    # from lme.
    probt1 <- pt(t1, df, lower.tail = T)
    probt2 <- pt(t2, df, lower.tail = F)
    ifelse(probt1 < alpha && probt2 < alpha, check_equiv <- "There is evidence in favour of a practically meaningless trend.", 
        check_equiv <- "There is no evidence in favour of practically meaninglesls trend.")
    res <- list(t1 = t1, t2 = t2, pValue.t1 = probt1, 
        pValue.t2 = probt2, decision = check_equiv, 
        ei = ei)
    class(res) <- "negTrend"
    return(res)
}

#' @rdname negTrend
#' @param x object of class \code{negTrend}
#' @export
print.negTrend <- function(x, ...) {
    cat("-------Two one-sided test for negligible trend------\n\n")
    cat("The equivalence interval magnitude was ", 
        x[[6]], "\n\n")
    cat("t1 = ", x[[1]], "p-value = ", x[[3]], "\n")
    cat("t2 = ", x[[2]], "p-value = ", x[[4]], "\n")
    cat(x[[5]])
} 
