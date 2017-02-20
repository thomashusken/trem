#' Fit a zero-truncated recurrent events model
#'
#' \code{fitrem} fits a zero-truncated recurrent events model on the data
#' provided and calculates the point and interval estimate of the population
#' size
#'
#' This is a specialised function, as it requires the data to be in a
#' three-dimensional format (persons x time x variables). The default
#' optimization procedure is a Newton-Raphson method, while numerical
#' optimization via \code{optim} is also possible
#'
#'
#' @param data An array with at least three dimensions
#' @param bstart A vector of starting values
#' @return A coefficient matrix, point and interval estimate of the population
#'   size, loglikelihood and number of iterations
#' @export


fitrem <- function(data, bstart, ...) {

    events <- data[, , 1]
    X <- data[, , -1, drop = F]
    p <- dim(X)[3]
    if (missing(bstart)) {
        bfirst <- bstart <- rep(0, p)
    } else if (length(bstart) < p) {
        warning("Starting values not provided for all parameters, using 0 as start value for parameters not specified")
        bfirst <- c(bstart, rep(0, p - length(bstart)))
    } else if (length(bstart) > p) {
        warning("Too many starting values provided, discarding excess values")
        bfirst <- bstart[1:p]
    } else bfirst <- bstart
    time <- dim(X)[2]
    n <- dim(X)[1]
    b <- matrix(bfirst, p, 1)
    iter <- 0
    tol <- 1

    while (tol > 1e-07) {
        u <- exp(tensor::tensor(X, b, alongA = 3, alongB = 1))
        dim(u) <- c(n, time)
        U <- c(rowSums(u))
        p0 <- exp(-U)
        pCap <- 1 - p0
        eU <- exp(U)
        logl <- sum(log(u[events == 1])) - sum(U) - sum(log(1 - p0))
        upCap <- (u * p0)/pCap
        EupCap <- (events - u - upCap)
        db <- tensor::tensor(X, EupCap, alongA = c(1, 2), alongB = c(1, 2))
        XuX <- tensor::tensor(X, c(u) * X, alongA = c(1, 2), alongB = c(1, 2))
        ddb1 <- tensor::tensor(X * c((u)/(eU - 1)), X, alongA = c(1, 2), alongB = c(1, 2))
        uX <- apply(c(u) * X, c(1, 3), sum)
        ddb2 <- colSums(tensor::tensor(X * c(u * eU)/(eU - 1)^2, uX, alongA = c(1), alongB = c(1)))
        ddb <- -XuX - ddb1 + ddb2
        infor <- solve(ddb)
        bnew <- b - infor %*% db
        tol <- sum(abs(bnew - b))
        b <- bnew
        iter <- iter + 1
        print(logl)
    }

    se <- sqrt(diag(-infor))
    tvalue <- b/se
    pvalue <- 2 * pnorm(-abs(tvalue))

    nb <- length(b)
    coef <- matrix(c(b, se, tvalue, pvalue), nrow = nb, ncol = 4)
    colnames(coef) <- c("B", "SE", "t-value", "P-value")
    rownames(coef) <- dimnames(data)[[3]][-1]

    Nhat <- sum(1/pCap)
    if (p != 1) {
        var1 <- colSums((uX * p0)/pCap^2)
    } else {
        var1 <- sum((uX * p0)/pCap^2)
    }
    var1 <- -(t(var1) %*% infor %*% var1)
    var2 <- sum(p0/pCap^2)
    Nse <- sqrt(var1 + var2)
    CI <- cbind(Nhat - 1.96*Nse, Nhat + 1.96*Nse)
    AIC <- -2 * logl + 2 * length(b)
    Rij <- cumsum(colSums(u))

    output <- list(coef = coef, logl = logl, AIC = AIC, iter = iter, Nhat = Nhat, CI = CI,
        Nse = Nse, infor = infor)
    class(output) <- "trem"
    output
}

print.trem <- function(x, ...){
  printCoefmat(x$coef, has.Pvalue = T)
  cat("\n", "Nhat (95%-CI) ", x$Nhat, " (" ,x$CI[1]," - ", x$CI[2], ")", sep = "")
  cat("\n\n", "Estimation converged in", x$iter, "iterations")
}

summary.trem <- function(x, ...){
  printCoefmat(x$coef, has.Pvalue = T)
  cat("\n", "Nhat (95%-CI) ", x$Nhat, " (" ,x$CI[1]," - ", x$CI[2], ")", sep = "")
  cat("\n\n", "Loglikelihood:", x$logl)
  cat("\n\n", "AIC:", x$AIC)
  cat("\n\n", "Estimation converged in", x$iter, "iterations")

}


