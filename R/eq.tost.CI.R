#' Independent samples equivalence test using confidence intervals
#' 
#' The two one-sided test for independent samples can be expressed using the confidence inclusion principle. If the confidence interval around a mean difference is within the equivalence bounds, then the mean difference is considered to be practically meaningless. 
#' 
#' @aliases eq.tost.CI
#' @param dat an N x 3 matrix or data.frame containing raw data used to compute the 
#'   correlation matrix between variables. The input may also be a 1 x 3 vector of correlations
#'   (r12, r13, and r23, respectively) and requires a sample size input (N)
#' @param ei equivalence interval
#' @param n sample size when dat input is a vector of correlations
#' @param alpha desired alpha level
#' 
#' @return returns a list containing the p-value, confidence interval, and statistical decision
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export eq.tost.CI
#' @examples
#' \dontrun{
#' #raw data
#' set.seed(1234)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' dat <- data.frame(dv=c(x,y),  group=c(rep('g1', length(x)), rep('g2', length(y)) ))
#' 
#' # Plot raw data distribution
#' 
#'  if (!require(ggplot2)) install.packages(ggplot2); library(ggplot2)
#' ggplot(dat,aes(x=dv)) + 
#'   geom_histogram(data=subset(dat,group == 'g1'),fill = 'red', alpha = 0.2) +
#'   geom_histogram(data=subset(dat,group == 'g2'),fill = 'blue', alpha = 0.2)
#'   
#' # Run test
#' ei <- 0.50
#' eq.tost.CI(x,y, ei=ei, alpha=.05)
#' equivInfo <- data.frame(name='mean difference', meanDiff=equivObj$meanDiff, lowCI=equivObj$lowCI, highCI=equivObj$hiCI)
#' 
#' p <- ggplot(equivInfo, aes(x=name, y=meanDiff))+
#'   geom_pointrange(aes(ymin=lowCI, ymax=highCI)) +  xlab('comparison') +
#'   geom_hline(yintercept=c(-equivint, equivint), linetype='dashed', color='blue') + 
#'   coord_flip()
#' p
#' eq.tost.CI(x,y, ei=0.50, alpha=.05)
#' }
#' 
eq.tost.CI <- function(x, y, ei, alpha = 0.05, na.rm=FALSE) {
  ifelse (any(is.na(y)) | any(is.na(x)), missing <- TRUE, missing <- FALSE    )
  if(missing & na.rm==TRUE){
      x <- na.omit(x)
      y <- na.omit(y)
  }
   if (missing & na.rm==FALSE){
     stop(print('There are missing values.'))
   }
    se <- sqrt(((((length(x) - 1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2))/(length(x) + 
        length(y) - 2)) * (1/length(x) + 1/length(y)))
    num1 <- (mean(x) - mean(y) - ei)
    num2 <- (mean(x) - mean(y) + ei)
    meanDiff <- mean(x) - mean(y)
    t1 <- (mean(x) - mean(y) - ei)/se
    t2 <- (mean(x) - mean(y) + ei)/se
    dft <- length(x) + length(y) - 2
    
    probt1 <- pt(t1, dft, lower.tail = T)
    probt2 <- pt(t2, dft, lower.tail = F)
    ifelse(probt1 <= alpha & probt2 <= alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
        decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
    
    # by two CIs (1-alpha) find the two CIs for each of the mean diffs
    # get critical values
    t1Crit <- qt(1 - alpha, dft, lower.tail = T)
    t2Crit <- qt(1 - alpha, dft, lower.tail = F)
    
    # by 1 CIs (1-2alpha)
    tCrit <- qt(1 - 2 * alpha, dft)
    lowCI <- meanDiff - tCrit * se
    hiCI <- meanDiff + tCrit * se
    ifelse(lowCI > -ei & hiCI < ei, ci.decis <- "CI bounds are within EI. Reject in favour of equivalence", 
        ci.decis <- "CI bounds are outside EI. Dont reject in favour of equivalence")
    
    ciBounds <- c(lowCI, hiCI)
    names(ciBounds) <- c("lowCI", "highCI")
    means <- c(mean(x), mean(y))
    names(means) <- c("Mean Grp 1", "Mean Grp 2")
    sds <- c(sd(x), sd(y))
    names(sds) <- c("SD Grp 1", "SD Grp 2")
    ei <- (c(ei))
    names(ei) <- c("equivalence interval")
    tstats <- c(t1, t2)
    dfs <- c(dft, dft)
    pvals <- c(probt1, probt2)
    names(tstats) <- c("t1", "t2")
    names(dfs) <- c("dft1", "dft2")
    names(pvals) <- c("p_t1", "p_t2")
    res <- list(means = means, meanDiff = meanDiff, sds = sds, ei = ei, 
        tstats = tstats, dfs = dfs, pvals = pvals, decis = decis, ciBounds = ciBounds, 
        ci.decis = ci.decis, se = se, twoAlphaTCrit = tCrit, lowCI = lowCI, 
        hiCI = hiCI)
    class(res) <- "eq.tost.CI"
    return(res)
} 
