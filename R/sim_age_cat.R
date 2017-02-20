sim_age_cat <- function(k = c(0.1, 0.3, 0.3, 0.2, 0.1), N, time, ...){
  AgeCat <- t(rmultinom(N, 1, k))
  Age3D  <- array(NA, dim = c(N, time, dim(AgeCat)[2] - 1))

  pdum <- dim(AgeCat)[2] - 1
  MoveCat <- rbinom(N, 1, 0.20)

  for(p in 1:pdum){
    Age3D[, , p] <- AgeCat[, p + 1]
  }

  for(i in 1:N){
    if(MoveCat[i] == 1){
      MoveDay <- sample(2:time, 1)
      if(any(Age3D[i, 1, 1:(pdum - 1)] > 0)){ #Middle dummy categories
        which.p <- which(Age3D[i, 1, ] == 1)
        Age3D[i, MoveDay:time, which.p] <- 0
        Age3D[i, MoveDay:time, which.p + 1] <- 1
      } else if(all(Age3D[i, 1, ] == 0)){ #Lowest dummy category
        Age3D[i, MoveDay:time, 1] <- 1
      } #No change to highest dummy category
    }
  }
  return(Age3D)
}
