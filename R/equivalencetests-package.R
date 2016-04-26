#' Package for Equivalence Tests
#'
#' The equivalencetests package provides statistical tools for equivalence testing in the null hypothesis testing framework. 
#' Available statistical tests include the following:   

#' ASSOCIATION-BASED  
#' [eq.corr] One-sample lack of correlation/association;     
#' [equiv_rs] Two independent samples equivalence of correlations (equivalence-based independent samples)    
#' [equiv_drs] Dependent samples equivalence of correlations (equivalence-based dependent samples)  
#' 
#' REPEATED-MEASURES / DEPENDENT SAMPLES   
#' [eq.hotT2] k-dependent samples equivalence of repeated/multiple measures (multivariate approach)    
#' [pw.std] k-dependent samples equivalence of repeated/multiple measures (pairwise comparisons approach, using standardized metric)   
#' [pw.unstd] k-dependent samples equivalence of repeated/multiple measures (pairwise comparisons approach, using unstand)    
#' [equivTrend] test for lack of trend   
#' 
#' GROUP-BASED / INDEPENDENT-SAMPLES  
#' [eq.tost] Two independent samples equivalence test of means  
#' [eq.tost.CI] Two independent samples equivalence test of means using the confidence-interval inclusion approach  
#' [eq.1way.ww] k independent samples equivalence tests of means   
#' [eq.normcomp] normative comparison tests for evaluating clinical significance   
#' 
#' @name equivalencetests-package
#' @aliases equivalencetests
#' @docType package
#' @title Equivalence tests
#' @author Rob Cribbie \email{cribbie@@yorku.ca}, Victoria Ng \email{ng.victoria.ky@@gmail.com}, Alyssa Counsell \email{counsela@@yorku.ca} and
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @import stats plyr mvtnorm ggplot2
#' @importFrom MBESS cor2cov
#' @keywords package
NULL

#' Description of AMDA data
#'
#' Abbreviated Master Data File Final (May 16, 2006)-imputed.
#'
#'
#' @name AMDA
#' @docType data
#' @author Rob Cribbie \email{cribbie@@yorku.ca}
#' @keywords data
NULL 
