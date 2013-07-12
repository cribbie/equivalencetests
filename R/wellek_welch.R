#' Wellek-Welch equivalence test
#' 
#' This function computes the Wellek-Welch equivalence test for k-groups.
#' 
#' @aliases wellek_welch
#' @param dv a numeric vector for the continous dependent variable
#' @param group grouping variable that can be coerced to a \code{factor}
#' @param eps numeric value defining the size of the equivalance interval
#' @param alpha the appropriate alpha level
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#' @export wellek_welch
#' @examples
#' \dontrun{
#' 
#' dv <- rnorm(100)
#' group <- rep(paste0('G',1:4), each=25)
#' wellek_welch(dv, group, eps = 0.25)
#' 
#' }
wellek_welch <- function(dv, group, eps, alpha = .05) {
    ret <- list()
    group <- factor(group)
    size <- table(group)
    ng <- length(size)    
    psisq <- oneway.test(dv ~ group,var.equal=TRUE)$statistic*((ng-1)/mean(size))
    crit_psisq <- ((ng-1)/mean(size)) * 
        qf(p=alpha,df1=ng-1,df2=sum(size)-ng,ncp=mean(size)*eps^2)
    ret$Wellek_F <- c(stat=as.numeric(psisq), df1=ng-1, df2=sum(size)-ng, p=crit_psisq)
    
    ## Compute the Wellek-Welch F test for Equivalence ##
    
    psisq <- oneway.test(dv ~ group)$statistic*((ng-1)/mean(size))  
    crit_psisq<-((ng-1)/mean(size)) * 
        qf(p=alpha,df1=ng-1,df2=oneway.test(dv ~ group)$parameter[2],ncp=mean(size)*eps^2)
    ret$Wellek_Welch <- c(stat=as.numeric(psisq), df1=ng-1, 
                          df2=oneway.test(dv ~ group)$parameter[2], p=crit_psisq)
    ret    
} 
