#' Data format for srplv
#' \code{data_to_srplv} converts the dataset to the srplv layout.
#'
#' @param data a data frame with three columns, recording the system ID, time and status.
#' @param sys.col the index containing the system ID.
#' @param time.col the index containing the event times and end-of-observation times.
#' @param status.col the index containing censoring or failure status.
#' @param n.comp the number of components.
#'
#' @return A list with the following components:
#'    \item{n}{the number of systems in the fleet.}
#'    \item{n.comp}{the number of components in the systems.}
#'    \item{r}{a vector containing the number of failures for each system in the fleet.}
#'    \item{time}{a list containing the failure time for each system.}
#'    \item{cens}{a vector that records the end-of-observation time for each system.}
#'
#' @export

data_to_srplv <- function(data, sys.col, time.col, status.col, n.comp){
  data <- data.frame(ID=data[,sys.col], t.fail=data[,time.col], status=data[,status.col])
  systms <- factor(unique(data$ID))
  n.sys <- length(systms)
  time_tot  <- vector("list",n.sys)
  r_tot <- rep(0,n.sys)
  for(i in 1:n.sys){
    time_tot[[i]] <- data$t.fail[data$ID==systms[i]]
    r_tot[i] <- length(time_tot[[i]])
  }
  time <- vector("list",n.sys)
  cens <- rep(NA,n.sys)
  r <- rep(0,n.sys)
  for(i in 1:n.sys){
    if(r_tot[i]!=1){
      time[[i]] <- time_tot[[i]][1:(r_tot[i]-1)]
      cens[i] <- time_tot[[i]][r_tot[i]]
    } else{
      time[[i]] <- NA
      cens[i] <- time_tot[[i]][r_tot[i]]
    }
    r[i] <- r_tot[i]-1
  }
  Data <- NULL
  Data$n  <- n.sys
  Data$n.comp <- n.comp
  Data$r     <- r
  Data$time   <- time
  Data$cens   <- cens
  return(Data)
}
