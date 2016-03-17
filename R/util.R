winvar <- function(x, tr = 0.2, na.rm = FALSE) {
    # Compute the gamma Winsorized variance for the data in the
    # vector x.  tr is the amount of Winsorization which defaults
    # to .2.
    if (na.rm) 
        x <- x[!is.na(x)]
    y <- sort(x)
    n <- length(x)
    ibot <- floor(tr * n) + 1
    itop <- length(x) - ibot + 1
    xbot <- y[ibot]
    xtop <- y[itop]
    y <- ifelse(y <= xbot, xbot, y)
    y <- ifelse(y >= xtop, xtop, y)
    winvar <- var(y)
    winvar
} 
