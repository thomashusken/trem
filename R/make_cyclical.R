make_cyclical <- function(data, cyc = 365){
  n             <- dim(data)[1]
  time          <- dim(data)[2]
  cost          <- matrix(cos(((2*pi)/cyc)*(1:time)), n, time, byrow = T)
  sint          <- matrix(sin(((2*pi)/cyc)*(1:time)), n, time, byrow = T)
  out           <- abind::abind(cost, sint, along = 3)
  dimnames(out)[[3]] <- list(paste0("Bcos",cyc), paste0("Bsin",cyc))
  out
}


