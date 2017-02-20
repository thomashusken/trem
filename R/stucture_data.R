structure_data <- function(data, dcols, xcolnames, otime, ...) {
    n <- dim(data)[1]
    events <- make_event_matrix(data, dcols, ...)
    covmatrix <- make_cov_array(data, xcolnames, otime, ...)
    data <- abind(events, covmatrix)
    output <- data
    dimnames(output)[[3]] <- c("C", paste0(dimnames(covmatrix)[[3]]))
    return(output)
}


make_event_matrix <- function(data, dcols, startt, sel, ...) {
    if (missing(startt)) {
        startt <- min(data[, dcols[1]])
    }

    if (class(data[, dcols[1]]) %in% c("POSIXct", "Date")) {
        DaysOfYear <- apply(data[, dcols], 2, function(x) ymd(x) - startt)
    } else if (class(data[, dcols[1]]) %in% c("numeric", "integer")) {
        DaysOfYear <- data[, dcols]
    }

    n <- dim(data)[1]

    if (missing(sel)) {
        otime <- max(ymd(data[, dcols[1]])) - min(ymd(data[, dcols[1]]))
        Events <- matrix(0, n, otime)
        for (i in 1:n) {
            Events[i, as.numeric(DaysOfYear[i, ])] <- 1
        }
    } else {
        if (class(sel) == "Date") {
            num <- (ymd(sel) + 1) - startt
            Events <- matrix(0, n, (max(num) + 1) - min(num))
            for (j in 1:dim(DaysOfYear)[2]) {
                pick <- apply(DaysOfYear, 2, function(x) which(x %in% num))[[j]]
                for (i in pick) {
                  Events[i, as.numeric(DaysOfYear[i, j]) - (min(num - 1))] <- 1
                }
            }
        } else if (class(sel) != "Date")
            stop("sel should be of class Date")
    }
    output <- Events
}

make_cov_array <- function(data, xcolnames, otime, ...) {
    xcols <- as.formula(paste("~", paste(xcolnames, collapse = "+")))
    covs <- model.matrix(xcols, data = data)
    n <- dim(covs)[1]
    k <- dim(covs)[2]
    covmatrix <- array(0, dim = c(n, k, otime))
    dimnames(CovMatrix) <- list(NULL, c(dimnames(covs)[[2]]), paste0("d", 1:otime))
    covmatrix[, , ] <- covs
    covmatrix <- aperm(covmatrix, perm = c(1, 3, 2))
}



