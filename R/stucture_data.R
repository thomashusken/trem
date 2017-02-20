#' Structure recurrent events data to 3d format
#'
#' \code{structure_data} restructures a wide-format matrix to a 3d array, which can be passed on to \code{fitrem} to fit a model
#'
#' @param data a wide-format matrix
#' @param dcols column names containing capture dates
#' @param xcolnames column names containing time-invariant covariates
#' @param otime the observation window in days
#' @return a 3D-array of dimension n x otime x \code{length(xcolnames)}
#' @export

structure_data <- function(data, dcols, xcolnames, otime, ...) {
    n <- dim(data)[1]
    events <- make_event_matrix(data, dcols, ...)
    covmatrix <- make_cov_array(data, xcolnames, otime, ...)
    data <- abind::abind(events, covmatrix)
    output <- data
    dimnames(output)[[3]] <- c("C", paste0(dimnames(covmatrix)[[3]]))
    return(output)
}


make_event_matrix <- function(data, dcols, startt, sel, ...) {
    if (missing(startt)) {
        startt <- min(data[, dcols[1]])
    }

    if (class(data[, dcols[1]]) %in% c("POSIXct", "Date")) {
        days_of_year <- apply(data[, dcols], 2, function(x) ymd(x) - startt)
    } else if (class(data[, dcols[1]]) %in% c("numeric", "integer")) {
        days_of_year <- data[, dcols]
    }

    n <- dim(data)[1]

    if (missing(sel)) {
        otime <- max(ymd(data[, dcols[1]])) - min(ymd(data[, dcols[1]]))
        events <- matrix(0, n, otime)
        for (i in 1:n) {
            events[i, as.numeric(days_of_year[i, ])] <- 1
        }
    } else {
        if (class(sel) == "Date") {
            num <- (ymd(sel) + 1) - startt
            events <- matrix(0, n, (max(num) + 1) - min(num))
            for (j in 1:dim(days_of_year)[2]) {
                pick <- apply(days_of_year, 2, function(x) which(x %in% num))[[j]]
                for (i in pick) {
                  events[i, as.numeric(days_of_year[i, j]) - (min(num - 1))] <- 1
                }
            }
        } else if (class(sel) != "Date")
            stop("sel should be of class Date")
    }
    output <- events
}

make_cov_array <- function(data, xcolnames, otime, ...) {
    xcols <- as.formula(paste("~", paste(xcolnames, collapse = "+")))
    covs <- model.matrix(xcols, data = data)
    n <- dim(covs)[1]
    k <- dim(covs)[2]
    covmatrix <- array(0, dim = c(n, k, otime))
    dimnames(covmatrix) <- list(NULL, c(dimnames(covs)[[2]]), paste0("d", 1:otime))
    covmatrix[, , ] <- covs
    covmatrix <- aperm(covmatrix, perm = c(1, 3, 2))
}



