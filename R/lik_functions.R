#' Function lik calculates log likelihood contribution of the fleet of systems.
#'
lik <- function(par,Data,delta,distribution) {
  out <- 0
  for (i in 1:Data$n) {
    out <- out +  eval(parse(text=paste('liksyst_', distribution, '(par=par,r=Data$r[i],time=Data$time[[i]],d=delta[[i]],
                             m=Data$n.comp,cens=Data$cens[i])', sep='')))
  }
  return(out)
}

#' Function liksyst calculates log likelihood contribution of a single system.
## -------------------- If the chosen distribution is Weibull distribution ------------ #
liksyst_weibull <- function(par,r,time,d,m,cens) {
  if (r != 0) {
    aux.time <- rep(NA,r)
    aux.time[1] <- time[1]
    if (r >= 2) {
      for (k in 2:r) {
        cond <- d[k] == d[1:(k-1)]
        aux.time[k] <- ifelse(any(cond), time[k] - time[max((1:(k-1))[cond[1:(k-1)]])],
                              time[k])
      }
    }
    comps.fails <- setdiff(1:m,setdiff(1:m,d))
    m.fails <- length(comps.fails)
    time.end <- rep(0,m.fails)
    for (a in 1:m.fails) {
      x <- comps.fails[a]==d
      time.end[a] <- time[max(which(x))]
    }
    out <- (sum(dweibull(aux.time,shape=par[1],scale=par[2],log=TRUE))
            + sum(pweibull(cens-time.end,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))
            +(m-m.fails)*pweibull(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))

  } else {
    out <- m*pweibull(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE)
  }
  return(out)
}

## -------------------- If the chosen distribution is gamma distribution ------------ #
liksyst_gamma <- function(par,r,time,d,m,cens) {
  if (r != 0) {
    aux.time <- rep(NA,r)
    aux.time[1] <- time[1]
    if (r >= 2) {
      for (k in 2:r) {
        cond <- d[k] == d[1:(k-1)]
        aux.time[k] <- ifelse(any(cond), time[k] - time[max((1:(k-1))[cond[1:(k-1)]])],
                              time[k])
      }
    }
    comps.fails <- setdiff(1:m,setdiff(1:m,d))
    m.fails <- length(comps.fails)
    time.end <- rep(0,m.fails)
    for (a in 1:m.fails) {
      x <- comps.fails[a]==d
      time.end[a] <- time[max(which(x))]
    }

    out <- (sum(dgamma(aux.time,shape=par[1],scale=par[2],log=TRUE))
            + sum(pgamma(cens-time.end,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))
            +(m-m.fails)*pgamma(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))

  } else {
    out <- m*pgamma(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE)
  }
  return(out)
}

## -------------------- If the chosen distribution is lognormal distribution ------------ #
liksyst_lnorm <- function(par,r,time,d,m,cens) {
  if (r != 0) {
    aux.time <- rep(NA,r)
    aux.time[1] <- time[1]
    if (r >= 2) {
      for (k in 2:r) {
        cond <- d[k] == d[1:(k-1)]
        aux.time[k] <- ifelse(any(cond), time[k] - time[max((1:(k-1))[cond[1:(k-1)]])],
                              time[k])
      }
    }
    comps.fails <- setdiff(1:m,setdiff(1:m,d))
    m.fails <- length(comps.fails)
    time.end <- rep(0,m.fails)
    for (a in 1:m.fails) {
      x <- comps.fails[a]==d
      time.end[a] <- time[max(which(x))]
    }
    out <- (sum(dlnorm(aux.time,meanlog=par[1],sdlog=par[2],log=TRUE))
            + sum(plnorm(cens-time.end,meanlog=par[1],sdlog=par[2],log=TRUE,lower.tail=FALSE))
            +(m-m.fails)*plnorm(cens,meanlog=par[1],sdlog=par[2],log=TRUE,lower.tail=FALSE))

  } else {
    out <- m*plnorm(cens,meanlog=par[1],sdlog=par[2],log=TRUE,lower.tail=FALSE)
  }
  return(out)
}

## -------------------- If the chosen distribution is llogis distribution ------------ #
liksyst_llogis <- function(par,r,time,d,m,cens) {
  if (r != 0) {
    aux.time <- rep(NA,r)
    aux.time[1] <- time[1]
    if (r >= 2) {
      for (k in 2:r) {
        cond <- d[k] == d[1:(k-1)]
        aux.time[k] <- ifelse(any(cond), time[k] - time[max((1:(k-1))[cond[1:(k-1)]])],
                              time[k])
      }
    }
    comps.fails <- setdiff(1:m,setdiff(1:m,d))
    m.fails <- length(comps.fails)
    time.end <- rep(0,m.fails)
    for (a in 1:m.fails) {
      x <- comps.fails[a]==d
      time.end[a] <- time[max(which(x))]
    }

    out <- (sum(flexsurv::dllogis(aux.time,shape=par[1],scale=par[2],log=TRUE))
            + sum(flexsurv::pllogis(cens-time.end,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))
            +(m-m.fails)*flexsurv::pllogis(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE))

  } else {
    out <- m*flexsurv::pllogis(cens,shape=par[1],scale=par[2],log=TRUE,lower.tail=FALSE)
  }
  return(out)
}

