#' Trend test
#' This function takes as input the estimates obtained from a linear mixed model. Specifically, some fixed effect slope for time or some other ordered variable, its standard error, and the approximated degrees of freedom are used. 
#' This is essentially a TOST; instead of the mean difference as the parameter of interest, the slope is. 
#' The equivalence interval is then specified with the lower and upper limit bounds for what is considered to be minimum average slope deemed practically important. 
#' @param  ei equivalence interval
#' @param slope the estimated fixed effect slope obtained from some model
#' @param se the standard error of the slope
#' @param df the degrees of freedom obtained from the fitted model  
#' 
#' equivTrend()
equivTrend <- function(ei, slope, se, df, n, alpha = 0.05) {
    t1 <- (slope - ei)/se
    t2 <- (slope + ei)/se
    # dft <- summary(mod)$tTable[2,3] #this df is from lme.
    probt1 <- pt(t1, df, lower.tail = T)
    probt2 <- pt(t2, df, lower.tail = F)
    ifelse(probt1 < alpha && probt2 < alpha, check_equiv <- "There is evidence in favour of a practically meaningless trend.", check_equiv <- 'There is no evidence in favour of practically meaninglesls trend.')
    res <- list(t1=t1, t2=t2, p.t1=probt1, p.t2=probt2, check_equiv)
    return(res)
} 
