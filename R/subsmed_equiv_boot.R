#' Computes the equivalence and bootstrap-based test of substantial mediation
#' 
#' This function computes the equivalence and bootstrap-based test of substantial mediation 
#' outlined in: Mara, C. & Cribbie, R. A. (submitted). An equivalence-based procedure for 
#' assessing substantial mediation. Psychological Methods.
#'
#' This function determines whether a single mediator accounts for a substantial proportion 
#' of the relationship between an independent variable (IV) and a dependent variable (DV), 
#' where quantifying 'a substantial proportion' depends on the width of the equivalence 
#' interval. For a mediator to be said to account for a substantial proportion
#' of the relationship between the IV and DV, the raw relationship between IV and DV 
#' should be very similar to the reproduced IV-DV relationship (i.e., IV->M*M->DV). 
#' The equivalence interval defines how similar these relationships
#' must be in order for the mediation to be substantial.
#' 
#' @aliases subsmed_equiv_boot
#' @param IV independent variable
#' @param DV dependent variable
#' @param m mediator
#' @param ei equivalence interval
#' @param standardize logical; indicate whether the analysis should be done on the raw data or on 
#'   standardized data (i.e., all variables having a mean of 0 and standard deviation of 1)
#' @param nboot number of bootstraps
#' @param alpha the appropriate alpha level
#' @author Constance Mara \email{constancemara@@gmail.com}, 
#'   Rob Cribbie \email{cribbie@@yorku.ca}, and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export subsmed_equiv_boot
#' @references 
#' 
#' Mara, C. & Cribbie, R. A. (submitted). An equivalence-based procedure for assessing substantial mediation.
#' \emph{Psychological Methods}.
#' 
#' @examples
#' \dontrun{
#' # equivalence interval of .1, a standardized solution, 1000 bootstrap samples, 
#' #  and the default alpha level (.05) 
#' m <- rnorm(500)
#' x <- m + rnorm(500) 
#' y <- m + rnorm(500) 
#' subsmed_equiv_boot(x,y,m,ei=.1,nboot=1000,standardize=TRUE) 
#' }
#' 
subsmed_equiv_boot <- function(IV, DV, m, ei, standardize = FALSE, 
    nboot = 10000, alpha = 0.05) {
    x <- IV
    y <- DV
    x <- x[is.na(x) == F & is.na(y) == F & is.na(m) == F]
    y <- y[is.na(x) == F & is.na(y) == F & is.na(m) == F]
    m <- m[is.na(x) == F & is.na(y) == F & is.na(m) == F]
    n <- length(x)
    if (standardize == TRUE) {
        x <- scale(x)
        y <- scale(y)
        m <- scale(m)
    }
    id <- 1:n
    newdat <- data.frame(x, y, m)
    semmod <- structure(c("x->m", "m->y", "m<->m", "y<->y", "x<->x", 
        "bmx", "bym", "erm", "ery", "erx", NA, NA, NA, NA, NA), 
        .Dim = c(5L, 3L), class = "semmod")
    moments <- sem::rawMoments(~(-1) + m + x + y)
    semout <- sem::sem(semmod, moments, n)
    covyx1 <- lm(y ~ x)$coeff[2]
    covyx2 <- semout$coeff[1] * semout$coeff[2]
    bootresults_eq <- matrix(0, nrow = nboot, ncol = 2)
    input <- matrix(0, nrow = 3, ncol = 1)
    results <- matrix(0, nrow = 2, ncol = 2)
    covs <- matrix(0, nrow = 2, ncol = 1)
    for (i in 1:nboot) {
        a <- sample(id, n, replace = T)
        newdatb <- newdat[a, ]
        xb <- newdatb[, 1]
        yb <- newdatb[, 2]
        mb <- newdatb[, 3]
        semmodb <- structure(c("xb->mb", "mb->yb", "mb<->mb", "yb<->yb", 
            "xb<->xb", "bmx", "bym", "erm", "ery", "erx", NA, NA, 
            NA, NA, NA), .Dim = c(5L, 3L), class = "semmod")
        momentsb <- sem::rawMoments(~(-1) + mb + xb + yb)
        semoutb <- sem::sem(semmodb, momentsb, n)
        covyx1b <- lm(yb ~ xb)$coeff[2]
        covyx2b <- semoutb$coeff[1] * semoutb$coeff[2]
        bootresults_eq[i, 1] <- covyx1b - covyx2b + ei
        bootresults_eq[i, 2] <- covyx1b - covyx2b - ei
    }
    ifelse(quantile(bootresults_eq[, 1], alpha) > 0 & quantile(bootresults_eq[, 
        2], 1 - alpha) < 0, decis <- "It is concluded that the mediator explains a substantial proportion of the variability in the relationship between the predictor and outcome", 
        decis <- "It is concluded that the mediator does not explain a substantial proportion of the variability in the relationship between the predictor and outcome")
    
    title1 <- "An Equivalence and Boostrap-based Test of Substantial Mediation"
    title2 <- "Input Parameters"
    input[1, 1] <- alpha
    input[2, 1] <- ei
    input[3, 1] <- nboot
    rownames(input) <- c("alpha", "Equiv Interval", "# of Boot Samples")
    colnames(input) <- c("Value")
    title3 <- "Important Relationships"
    covs[1, 1] <- covyx1
    covs[2, 1] <- covyx2
    rownames(covs) <- c("raw X-Y relationship", "reproduced X-Y relationship, (X->M)*(M->Y)")
    colnames(covs) <- c("Value")
    title4 <- "Test Statistic Confidence Intervals"
    results[1, 1] <- quantile(bootresults_eq[, 1], alpha)
    results[1, 2] <- quantile(bootresults_eq[, 1], 1 - alpha)
    results[2, 1] <- quantile(bootresults_eq[, 2], alpha)
    results[2, 2] <- quantile(bootresults_eq[, 2], 1 - alpha)
    rownames(results) <- c("Positive Equiv Interval", "Negative Equiv Interval")
    colnames(results) <- c("alpha%CI", "1-alpha%CI")
    title5 <- "Decision"
    out <- list(title1, title2, input, title3, covs, title4, results, 
        title5, decis)
    class(out) <- "subsmed_equiv_boot"
    out
}

#' @S3method print subsmed_equiv_boot
#' @rdname subsmed_equiv_boot
#' @method print subsmed_equiv_boot
#' @param x object of class \code{subsmed_equiv_boot}
#' @param ... additional arguments
print.subsmed_equiv_boot <- function(x, ...) {
    cat(x[[1]], "\n\n")
    cat(x[[2]], "\n")
    print(x[[3]])
    cat("\n", x[[4]], "\n")
    print(x[[5]])
    cat("\n", x[[6]], "\n")
    print(x[[7]])
    cat("\n", x[[9]], "\n")
} 
