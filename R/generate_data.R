#' Data generation of Superposed Renewal Processes
#'
#' \code{generate_data} generates srplv dataset.
#'
#' @param n.sys  The number of systems.
#' @param n.comp  The number of sockets for each system.
#' @param m.t  The expected failure time of components.
#' @param v.t The variance of failure time of components .
#' @param m.c The expected failure time of censor (end-of-observation time).
#' @param v.c  The variance of failure time of censor (end-of-observation time).
#' @param distribution The chose distribution: weibull, gamma, lnorm, llogis.
#' @param s.seed  The seed for data generation.
#'
#' @return A list with the following components:
#'    \item{n}{the number of systems in the fleet.}
#'    \item{n.comp}{the number of sockets in each system.}
#'    \item{r}{a vector with a length n, recording the number of failures for each system.}
#'    \item{time}{a list that each vector records the failure times of each system.}
#'    \item{cens}{a vector with a length n, recording the end-of-observation time for each system.}
#'    \item{par.t}{a vector with parameters of components failure times distribution.}
#'    \item{par.c}{a vector with parameters of censor distribution.}
#'    \item{delta.true}{a list that each vector records the socket index failure indicator of each system.}
#'
#' @export
#'
generate_data <- function(n.sys,n.comp,m.t,v.t,m.c,v.c,distribution,s.seed){
  out <- try(eval(parse(text=paste('ger_', distribution, '(n.sys=n.sys,n.comp=n.comp,m.t=m.t,v.t=v.t,m.c=m.c,v.c=v.c,s.seed=s.seed)', sep=''))),silent=TRUE)
  if(class(out)=="try-error"){
    print("Error - please select valid arguments.")
  } else {
    return(out)
  }
}

# --------------------------------------------------------------------------------------------- #
# -------------------- Functions with mean (mu) and variance (s2) as input and ---------------- #
# ----------------------- to return distibution parameters ------------------------------------ #

## -------------------- Function to return shape and scale of Weibull distribution ------------ #
fun_weibull <- function(mu,s2){
  f.sh <- function(x) {log((s2+mu^2)/mu^2) - lgamma(1+2/x) + 2*lgamma(1+1/x)}
  sh <- uniroot(f.sh,interval=c(10^-200, 10^10))$root
  scl <- mu/gamma(1+(1/sh))
  return(c(sh,scl))
}
## -------------------- Generate data from Weibull distribution ------------------------------- #
ger_weibull <- function(n.sys,n.comp,m.t,v.t,m.c,v.c,s.seed){
  delta <- vector("list",n.sys)
  time  <- vector("list",n.sys)
  r <- rep(0,n.sys)
  par.c <- fun_weibull(mu=m.c,s2=v.c)
  par.t <- fun_weibull(mu=m.t,s2=v.t)
  set.seed(s.seed)
  cens <- rweibull(n.sys,shape=par.c[1],scale=par.c[2])
  for (i in 1:n.sys) {
    j <- 1
    fail  <- rweibull(n.comp,shape=par.t[1],scale=par.t[2])
    while(min(fail) < cens[i]) {
      delta[[i]] <- c(delta[[i]],order(fail)[1])
      aux.f <- fail
      fail[delta[[i]][j]] <- fail[delta[[i]][j]]+rweibull(1,shape=par.t[1],scale=par.t[2])
      time[[i]] <- c(time[[i]],min(aux.f))
      j <- j+1
    }
    r[i] <- length(delta[[i]])
  }
  tau <- vector("list",n.sys)
  for (i in 1:n.sys){
    if(!is.null(time[[i]])){
      tau[[i]][1] <- time[[i]][1]
      if(length(time[[i]])>1){
        for (k in 2:length(time[[i]])){
          tau[[i]][k] <- time[[i]][k] - time[[i]][k-1]
        }
      }
    }
  }
  Data <- NULL
  Data$n  <- n.sys
  Data$n.comp <- n.comp
  Data$r     <- r
  Data$time   <- time
  Data$cens   <- cens
  Data$par.t <- par.t
  Data$par.c <- par.c
  Data$delta.true <- delta
  return(Data)
}

## -------------------- Function to return shape and scale of gamma distribution ------------ #
fun_gamma <- function(mu,s2){
  sh <- mu^2/s2
  scl  <- s2/mu
  return(c(sh,scl))
}
## -------------------- Generate data from gamma distribution ------------------------------- #
ger_gamma <- function(n.sys,n.comp,m.t,v.t,m.c,v.c,s.seed){
  delta <- vector("list",n.sys)
  time  <- vector("list",n.sys)
  r <- rep(0,n.sys)
  par.c <- fun_gamma(mu=m.c,s2=v.c)
  par.t <- fun_gamma(mu=m.t,s2=v.t)
  set.seed(s.seed)
  cens <- rgamma(n.sys,shape=par.c[1],scale=par.c[2])
  for (i in 1:n.sys) {
    j <- 1
    fail  <- rgamma(n.comp,shape=par.t[1],scale=par.t[2])
    while(min(fail) < cens[i]) {
      delta[[i]] <- c(delta[[i]],order(fail)[1])
      aux.f <- fail
      fail[delta[[i]][j]] <- fail[delta[[i]][j]]+rgamma(1,shape=par.t[1],scale=par.t[2])
      time[[i]] <- c(time[[i]],min(aux.f))
      j <- j+1
    }
    r[i] <- length(delta[[i]])
  }
  tau <- vector("list",n.sys)
  for (i in 1:n.sys){
    if(!is.null(time[[i]])){
      tau[[i]][1] <- time[[i]][1]
      if(length(time[[i]])>1){
        for (k in 2:length(time[[i]])){
          tau[[i]][k] <- time[[i]][k] - time[[i]][k-1]
        }
      }
    }
  }
  Data <- NULL
  Data$n  <- n.sys
  Data$n.comp <- n.comp
  Data$r     <- r
  Data$time   <- time
  Data$cens   <- cens
  Data$par.t <- par.t
  Data$par.c <- par.c
  Data$delta.true <- delta
  return(Data)
}

## -------------------- Function to return meanlog and sdlog of lognormal distribution ------------ #
fun_lnorm <- function(mu,s2){
  ml  <- log(mu^2/sqrt(mu^2+s2))
  sl <- sqrt(log((mu^2+s2)/mu^2))
  return(c(ml,sl))
}
## -------------------- Generate data from lognormal distribution ------------------------------- #
ger_lnorm <- function(n.sys,n.comp,m.t,v.t,m.c,v.c,s.seed){
  delta <- vector("list",n.sys)
  time  <- vector("list",n.sys)
  r <- rep(0,n.sys)
  par.c <- fun_lnorm(mu=m.c,s2=v.c)
  par.t <- fun_lnorm(mu=m.t,s2=v.t)
  set.seed(s.seed)
  cens <- rlnorm(n.sys,meanlog=par.c[1],sdlog=par.c[2])
  for (i in 1:n.sys) {
    j <- 1
    fail  <- rlnorm(n.comp,meanlog=par.t[1],sdlog=par.t[2])
    while(min(fail) < cens[i]) {
      delta[[i]] <- c(delta[[i]],order(fail)[1])
      aux.f <- fail
      fail[delta[[i]][j]] <- fail[delta[[i]][j]]+rlnorm(1,meanlog=par.t[1],sdlog=par.t[2])
      time[[i]] <- c(time[[i]],min(aux.f))
      j <- j+1
    }
    r[i] <- length(delta[[i]])
  }
  tau <- vector("list",n.sys)
  for (i in 1:n.sys){
    if(!is.null(time[[i]])){
      tau[[i]][1] <- time[[i]][1]
      if(length(time[[i]])>1){
        for (k in 2:length(time[[i]])){
          tau[[i]][k] <- time[[i]][k] - time[[i]][k-1]
        }
      }
    }
  }
  Data <- NULL
  Data$n  <- n.sys
  Data$n.comp <- n.comp
  Data$r     <- r
  Data$time   <- time
  Data$cens   <- cens
  Data$par.t <- par.t
  Data$par.c <- par.c
  Data$delta.true <- delta
  return(Data)
}

## ------------- Function to return meanlog and sdlog of inverse gaussian distribution ---------- #
fun_invgauss <- function(mu,s2){
  m  <- mu
  disp <- s2/(mu^3)
  return(c(m,disp))
}
## -------------------- Generate data from inverse gaussian distribution ------------------------- #
ger_invgauss <- function(n.sys,n.comp,m.t,v.t,m.c,v.c,s.seed){
  delta <- vector("list",n.sys)
  time  <- vector("list",n.sys)
  r <- rep(0,n.sys)
  par.c <- fun_invgauss(mu=m.c,s2=v.c)
  par.t <- fun_invgauss(mu=m.t,s2=v.t)
  set.seed(s.seed)
  cens <- statmod::rinvgauss(n.sys,mean=par.c[1],dispersion=par.c[2])
  for (i in 1:n.sys) {
    j <- 1
    fail  <- statmod::rinvgauss(n.comp,mean=par.t[1],dispersion=par.t[2])
    while(min(fail) < cens[i]) {
      delta[[i]] <- c(delta[[i]],order(fail)[1])
      aux.f <- fail
      fail[delta[[i]][j]] <- fail[delta[[i]][j]]+statmod::rinvgauss(1,mean=par.t[1],dispersion=par.t[2])
      time[[i]] <- c(time[[i]],min(aux.f))
      j <- j+1
    }
    r[i] <- length(delta[[i]])
  }
  tau <- vector("list",n.sys)
  for (i in 1:n.sys){
    if(!is.null(time[[i]])){
      tau[[i]][1] <- time[[i]][1]
      if(length(time[[i]])>1){
        for (k in 2:length(time[[i]])){
          tau[[i]][k] <- time[[i]][k] - time[[i]][k-1]
        }
      }
    }
  }
  Data <- NULL
  Data$n  <- n.sys
  Data$n.comp <- n.comp
  Data$r     <- r
  Data$time   <- time
  Data$cens   <- cens
  Data$par.t <- par.t
  Data$par.c <- par.c
  Data$delta.true <- delta
  return(Data)
}


