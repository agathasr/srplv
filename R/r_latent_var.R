#' Latent variables generation
#'
#' \code{rdelta} generates latent variables for all systems in the fleet.
#'
#' @param par  The vector of parameters.
#' @param Data  Data in the srplv format as \code{\link{data_to_srplv}}.
#' @param distribution The chose distribution: weibull, gamma, lnorm, llogis.
#'
#' @return A list with the latent variables for all systems in the fleet.
#' @export
#'
rdelta <- function(par,Data,distribution) {
  if(distribution=='lnorm'){
    par[2] <- exp(par[2])
  } else{
    par <- exp(par)
  }
  delta <- vector("list",Data$n)
  for (i in 1:Data$n) {
    if(Data$r[i] != 0) {
      delta[[i]] <- eval(parse(text=paste('gera_d_', distribution, '(par=par,m=Data$n.comp,r=Data$r[i],time=Data$time[[i]])', sep='')))
    }
    else{
      delta[[i]] <- 0
    }
  }
  return(delta)
}


# --------------------------------------------------------------------------------------------- #
# -------------------- Functions to generate the latent variables ----------------------------- #
# --------------------------------------------------------------------------------------------- #
## -------------------- If the chosen distribution is the weibull distribution ------------ #
gera_d_weibull <- function(par,m,r,time){
  aux.time <- rep(0,m)
  delta <- rep(NA,r)
  for (k in 1:r) {
    w <- dweibull(time[k]-aux.time,shape=par[1],scale=par[2])
    w <- w/sum(w)
    delta[k] <- (1:m)[rmultinom(1,1,w) == 1]
    aux.time[delta[k]] <- time[k]
  }
  return(delta)
}

## -------------------- If the chosen distribution is the gamma distribution ------------ #
gera_d_gamma <- function(par,m,r,time){
  aux.time <- rep(0,m)
  delta <- rep(NA,r)
  for (k in 1:r) {
    w <- dgamma(time[k]-aux.time,shape=par[1],scale=par[2])
    w <- w/sum(w)
    delta[k] <- (1:m)[rmultinom(1,1,w) == 1]
    aux.time[delta[k]] <- time[k]
  }
  return(delta)
}

## -------------------- If the chosen distribution is the lognormal distribution ------------ #
gera_d_lnorm <- function(par,m,r,time){
  aux.time <- rep(0,m)
  delta <- rep(NA,r)
  for (k in 1:r) {
    w <- dlnorm(time[k]-aux.time,meanlog=par[1],sdlog=par[2])
    w <- w/sum(w)
    delta[k] <- (1:m)[rmultinom(1,1,w) == 1]
    aux.time[delta[k]] <- time[k]
  }
  return(delta)
}

## -------------------- If the chosen distribution is the llogis distribution ------------ #
gera_d_llogis <- function(par,m,r,time){
  aux.time <- rep(0,m)
  delta <- rep(NA,r)
  for (k in 1:r) {
    w <- flexsurv::dllogis(time[k]-aux.time,shape=par[1],scale=par[2])
    w <- w/sum(w)
    delta[k] <- (1:m)[rmultinom(1,1,w) == 1]
    aux.time[delta[k]] <- time[k]
  }
  return(delta)
}

