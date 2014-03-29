#' Normative comparison tests for evaluating clinical significance
#'
#' The following function conducts normative comparison tests for evaluating clinical 
#' significance, using both the Cribbie and Arpin-Cribbie (2009) and Kendall et al. (1999) 
#' methods. Normative comparison tests compare a treated sample to a normal comparison 
#' sample to determine if they are statistically equivalent.
#' 
#' @aliases normcomp
#' @param premean pretest mean
#' @param presd pretest standard deviation
#' @param pren pretest sample size
#' @param postmean posttest mean
#' @param postsd post test standard deviation
#' @param postn post test sample size
#' @param normmean normal comparison mean
#' @param normsd normal comparison standard deviation
#' @param normn normal comparison sample size
#' @param alpha desired alpha level
#' @param ... additional arguments to be passed
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export normcomp
#' @references 
#' 
#' Cribbie, R. A. & Arpin-Cribbie, C. A. (2009). Evaluating clinical significance through equivalence testing:
#' Extending the normative comparisons approach. \emph{Psychotherapy Research, 19}, 677-686.
#' 
#' Kendall, P. C., Marrs-Garcia, A., Nath, S. R., & Sheldrick, R. C. (1999). Normative comparisons for the evaluation
#' of clinical significance. \emph{Journal of Consulting and Clinical Psychology, 67}, 285-299.
#' 
#' @examples
#' \dontrun{
#' #In this example the pretest mean was 12, pretest sd was 4, pretest sample size was 100, posttest mean was 8,
#' #posttest sd was 5, posttest sample size was 95, normal comparison mean was 7, normal comparison sd was 3,
#' #normal comparison sample size was 500. 
#' #The desired alpha level is .05, which is the default, so no alpha level is specified
#' normcomp(12,4,100,8,5,95,7,3,500)
#' }
normcomp <- function(premean, presd, pren, postmean, postsd, postn, normmean, normsd, normn, alpha = 0.05, 
    ...) {
    equivint <- c(0.5 * (normsd), normsd, 1.5 * (normsd))
    t1 <- c(0, 0, 0)
    t2 <- c(0, 0, 0)
    probt1 <- c(0, 0, 0)
    probt2 <- c(0, 0, 0)
    decis <- c(0, 0, 0)
    preresults <- data.frame(matrix(0, nrow = 1, ncol = 3, dimnames = list(c("Value"), c("t-test", "df", 
        "p-value"))))
    results <- matrix(0, nrow = 6, ncol = 3, dimnames = list(c("t-test1", "t-test2", "df1", "df2", "pvalue1", 
        "pvalue2"), c("EI=.5*sd", "EI=sd", "EI=1.5*sd")))
    results_kendall <- matrix(0, nrow = 6, ncol = 1, dimnames = list(c("t-test1", "t-test2", "df1", "df2", 
        "pvalue1", "pvalue2"), c("EI=sd")))
    pret <- (premean - normmean)/sqrt((presd^2/pren) + (normsd^2/normn))
    df_pret <- (((presd^2/pren) + (normsd^2/normn))^2)/((presd^4/(pren^2 * (pren - 1))) + (normsd^4/(normn^2 * 
        (normn - 1))))
    pval_pret <- pt(abs(pret), df_pret, lower.tail = F)
    ifelse(pval_pret < alpha, decis_pret <- "The pretest mean differs from the normative mean and the following normative comparisons are meaningful", 
        decis_pret <- "The pretest was not found to differ from the normative mean and therefore the following normative comparisons are not meaningful")
    for (i in 1:length(equivint)) {
        t1[i] <- (postmean - normmean - equivint[i])/sqrt((postsd^2/postn) + (normsd^2/normn))
        t2[i] <- (postmean - normmean + equivint[i])/sqrt((postsd^2/postn) + (normsd^2/normn))
        dft <- (((postsd^2/postn) + (normsd^2/normn))^2)/((postsd^4/(postn^2 * (postn - 1))) + (normsd^4/(normn^2 * 
            (normn - 1))))
        probt1[i] <- pt(t1[i], dft, lower.tail = T)
        probt2[i] <- pt(t2[i], dft, lower.tail = F)
        ifelse(probt1[i] < alpha & probt2[i] < alpha, decis[i] <- "The treated and normative means are declared equivalent at this Equivalence Interval", 
            decis[i] <- "The treated and normative means cannot be declared equivalent at this Equivalence Interval")
    }
    tk1 <- (postmean - normmean - equivint[2])/sqrt(((((postn - 1) * postsd^2) + ((normn - 1) * normsd^2))/(postn + 
        normn - 2)) * (1/postn + 1/normn))
    tk2 <- (postmean - normmean + equivint[2])/sqrt(((((postn - 1) * postsd^2) + ((normn - 1) * normsd^2))/(postn + 
        normn - 2)) * (1/postn + 1/normn))
    dfk1 <- postn + normn - 2
    dfk2 <- postn + normn - 2
    pvalk1 <- pt(tk1, dfk1, lower.tail = T)
    pvalk2 <- pt(tk2, dfk2, lower.tail = F)
    ifelse(pvalk1 < alpha & pvalk2 < alpha, decisk <- "The treated and normative means are declared equivalent with the Kendall et al. method using an Equivalence Interval of one SD of the normative group", 
        decisk <- "The treated and normative means cannot be declared equivalent with the Kendall method using an Equivalence Interval of one SD of the normative group")
    title1 <- "Normative Comparison Tests for Assessing Clinical Significance using the Cribbie and Arpin-Cribbie (2009) and Kendall et al. (1999) Procedures"
    title2 <- "Comparison of Pretest Results to Normative Data"
    preresults[1, 1] <- pret
    preresults[1, 2] <- df_pret
    preresults[1, 3] <- pval_pret
    # preresults[1, 4] <- decis_pret
    title3 <- "Test Statistics for the Cribbie & Arpin-Cribbie Method at Equivalence Intervals (EI) of .5, 1 and 1.5 times the SD of the Normal Comparison Group"
    for (j in 1:length(equivint)) {
        results[1, j] <- t1[j]
        results[2, j] <- t2[j]
        results[3, j] <- dft
        results[4, j] <- dft
        results[5, j] <- probt1[j]
        results[6, j] <- probt2[j]
    }
    title4 <- "Conclusions from the Cribbie & Arpin-Cribbie Method"
    title5 <- "Results for the Kendall et al. Normative Comparisons Method with the Equivalence Interval set to the SD of the Normal Comparison Group"
    title6 <- "Note that the results for the Kendall et. al method will not be accurate if the standard deviations of the treated and normal comparison groups are not equal"
    results_kendall[1, 1] <- tk1
    results_kendall[2, 1] <- tk2
    results_kendall[3, 1] <- dfk1
    results_kendall[4, 1] <- dfk2
    results_kendall[5, 1] <- pvalk1
    results_kendall[6, 1] <- pvalk2
    # results_kendall[7, 1] <- decisk
    decism <- as.matrix(decis)
    rownames(decism) <- c("EI=.5*sd", "EI=sd", "EI=1.5*sd")
    rownames(results_kendall) <- c("t-test1", "t-test2", "df1", "df2", "pvalue1", "pvalue2")
    out <- list(title1, title2, preresults, decis_pret, title3, results, title4, decism, title5, title6, 
        results_kendall, decisk)
    class(out) <- "normcomp"
    out
}

#' @S3method print normcomp
#' @rdname normcomp
#' @method print normcomp
#' @param x object of class \code{normcomp}
print.normcomp <- function(x, ...) {
    cat(x[[1]], "\n")
    cat("-----------------------------------------------------\n\n")
    cat(x[[2]], "\n")
    print(x[[3]])
    cat(x[[4]], "\n\n")
    cat(x[[5]], "\n")
    print(round(x[[6]], 4))
    cat(x[[7]], "\n")
    colnames(x[[8]]) <- ""
    print(x[[8]])
    cat("\n\n", x[[9]], "\n")
    cat(x[[10]], "\n")
    print(round(x[[11]], 4))
    cat(x[[12]])
} 
