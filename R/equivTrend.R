#' Trend test
#' This function takes as input the estimates obtained from a linear mixed model (fixed slope, random intercept model). 
#' This is essentially a TOST; instead of the mean difference as the parameter of interest, the slope is. 
#' @param  ei equivalence interval
#' @param slope the estimated fixed effect slope obtained from some model
#' @param se the standard error of the slope
#' @param df the degrees of freedom obtained from the fitted model  
#' equivTrend()
equivTrend <- function(ei, slope, se, df, n, alpha = 0.05) {
    t1 <- (slope - ei)/se
    t2 <- (slope + ei)/se
    # dft <- summary(mod)$tTable[2,3] #this df is
    # from lme.
    probt1 <- pt(t1, df, lower.tail = T)
    probt2 <- pt(t2, df, lower.tail = F)
    ifelse(probt1 < alpha && probt2 < alpha, check_equiv <- 1, 
        check_equiv <- 0)
    return(check_equiv)
}
 
