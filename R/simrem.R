#'Simulate zero-truncated recurrent events data
#'
#'\code{simrem} simulates zero-truncated events data and returns a 3D-array
#'
#'
#'@param N the total (non-truncated) population size
#'@param time the observation window
#'@param b0 the intercept, default value is -7
#'@param bpers coefficients for time-invariant covariates. The data is simulated as ~ B(0.5)
#'@param bage coefficients for time-varying dummy variables.
#'@param cyc logical, whether cyclical predictors are added to the data. Default is \code{FALSE}
#'@param alpha the amplitude of the cyclical effect. Only used when \code{cyc} is \code{TRUE}
#'@param theta the phase of the cyclical effect. Only used when \code{cyc} is \code{TRUE}
#'@param period the period of the cyclical effect, defaults to \code{time}
#'@param int specify whether an interaction between the cyclical effect and the first variable of \code{bpers} should be included
#'@param aint the amplitude of the cyclical interaction effect
#'@param tint the phase of the cyclical interaction effect
#'@return a 3D-array of dimensions n x time x variables. Note that n will vary over each simulation due to sampling fluctuation by drawing from
#'the total population N
#'@export



simrem <- function(N=5000, time=365, b0=-7, bpers=NULL, bage = NULL, cyc = F, alpha=NULL, theta=NULL, period=time, int = F, aint=NULL, tint = NULL, ...){
  if(!is.numeric(bage) & !is.null(bage)) stop("Invalid input for age coefficient")
  b        <- c(b0,bpers,bage)
  lbpers   <- length(bpers[bpers!=0])
  lbage    <- length(bage)
  if(cyc == T){
    bcos     <- alpha*cos(theta)
    bsin     <- alpha*sin(theta)
    bcyc     <- c(bcos, bsin)
    b        <- c(b, bcyc)
    lbcyc    <- length(bcyc)
    if(int == T){
      bcosint  <- aint*cos(tint) - alpha*cos(theta)
      bsinint  <- aint*sin(tint) - alpha*sin(theta)
      bcycint  <- c(bcosint, bsinint)
      b        <- c(b, bcycint)
      lbint    <- length(bcycint)
    }
  }


  nb       <- length(b)
  X        <- array(NA, dim=c(N, time, 1))
  X[, , 1] <- 1

  if(any(bpers!= 0)){
    pers <- aperm(array(rbinom(N*lbpers, 1, 0.5),  dim = c(N, lbpers, time)), c(1,3,2))
    X    <- abind::abind(X, pers)
  }

  if(!is.null(bage)){
    age  <- sim_age_cat(k = rep(0.2, lbage + 1), N = N, time = time)
    X    <- abind::abind(X, age)
  }

  if(cyc == T & !is.null(alpha)){
    for(p in 1:length(period)){
      Periodic <- make_cyclical(X, cyc = period[p])
      X             <- abind::abind(X, Periodic, along = 3)
    }
  }

  if(int == T & !is.null(aint)){
    dimsin  <- dim(X)[3]
    dimcos  <- dimsin - 1
    dimpers <- 1 + lbpers
    costint <- X[, , dimcos] * X[, , dimpers]
    sintint <- X[, , dimsin] * X[, , dimpers]
    X       <- abind::abind(X, costint, sintint)
  }


  u      <- exp(tensor::tensor(X, b, alongA = 3, alongB = 1))
  E      <- matrix(0, N, time)

  for(i in 1:N){
    for(j in 1:time){
      if(runif(1) < u[i,j]) E[i,j] <- 1
    }
  }


  NE       <- matrix(rowSums(E),ncol=1)
  X        <- X[NE>0,,]
  E        <- E[NE>0,]
  NE       <- NE[NE>0,]
  output   <- abind::abind(E,X,along=3)

  if(is.null(bage)){
    if(cyc == T & int == F){
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)), paste0(c("Bcos", "Bsin"), rep(period, each = 2)))
    } else if(cyc == T & int == T){
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)), paste0(c("Bcos", "Bsin"), rep(period, each = 2)),
                                 paste0(c("Bcos", "Bsin"), rep(period, each = 2), "x"))
    } else {
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)))
    }
  } else if(!is.null(bage)){
    if(cyc == T & int == F){
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)), paste0("bage",1:(lbage)), paste0(c("Bcos", "Bsin"), rep(period, each = 2)))
    } else if(cyc == T & int == T){
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)), paste0("bage",1:(lbage)), paste0(c("Bcos", "Bsin"), rep(period, each = 2)),
                                 paste0(c("Bcos", "Bsin"), rep(period, each = 2), "x"))
    } else {
      dimnames(output)[[3]] <- c("Y", paste0("b",0:(lbpers)), paste0("bage",1:(lbage)))
    }
  }
  return(output)
}

