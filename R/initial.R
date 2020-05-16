#' Initial values for the parameters
#'
#' \code{initials}  obtains initial values for the parameters.
#'
#' @param Data  Data in the srplv format as \code{\link{data_to_srplv}}.
#' @param distribution The chose distribution: weibull, gamma, lnorm, llogis.
#' @return A vector with initial parameters values.
#' @export
#'

initials <- function(Data,distribution){
  out <- eval(parse(text=paste('initials.', distribution, '(Data)',sep='')))
  return(out)
}

###### Data manipulation for flexsurvreg function used to generate initial values for parameters
initials_data <- function(Data){
  ind.greater <-which(Data$r>Data$n.comp)
  ind.lower <- setdiff(1:Data$n,ind.greater)
  data.aux <- NULL
  aux.list <- vector("list",length(ind.lower))
  if(is.null(ind.lower)==FALSE){
    for (i in ind.lower){
      if(Data$r[i]==0){
        data.aux <- rbind(data.aux,cbind(rep(Data$cens[i],Data$n.comp),rep(0,Data$n.comp)))
      } else{
        aux.list[[i]] <- Data$cens[i] - Data$time[[i]]
        data.aux <- rbind(data.aux, cbind(unlist(Data$time[i]),rep(1,length(Data$time[i]))),cbind(unlist(aux.list[i]),rep(0,length(aux.list[i]))),cbind(rep(Data$cens[i],Data$n.comp-Data$r[i]),rep(0,Data$n.comp-Data$r[i])))
      }
    }
  }
  if(is.null(ind.greater)==FALSE){
    for( i in ind.greater){
      delta <- rep(0,Data$r[i])
      for (j in 1:Data$r[i]){
        delta[j] <- (1:Data$n.comp)[rmultinom(1,1,rep(1/Data$n.comp,Data$n.comp)) == 1]
      }
      time <- unlist(Data$time[i])
      aux.time <- rep(NA,Data$r[i])
      aux.time[1] <- time[1]
      if (Data$r[i] >= 2) {
        for (k in 2:Data$r[i]) {
          cond <- delta[k] == delta[1:(k-1)]
          aux.time[k] <- ifelse(any(cond), time[k] - time[max((1:(k-1))[cond[1:(k-1)]])],
                                time[k])
        }
      }
      comps.fails <- setdiff(1:Data$n.comp,setdiff(1:Data$n.comp,delta))
      m.fails <- length(comps.fails)
      time.end <- rep(0,m.fails)
      for (a in 1:m.fails) {
        x <- comps.fails[a]==delta
        time.end[a] <- time[max(which(x))]
      }

      data.aux <- rbind(data.aux,cbind(aux.time,rep(1,length(aux.time))),cbind(Data$cens[i]-time.end,rep(0,length(time.end))),
                        cbind(rep(Data$cens[i],Data$n.comp-m.fails),rep(0,Data$n.comp-m.fails)))
    }
  }
  Data2 <- data.frame(t=data.aux[,1],status=data.aux[,2])
  return(Data2)
}

#### Initial valures for the Weibull parameters
initials.weibull <- function(Data){
  Data2 <- initials_data(Data)
  survfit <- suppressWarnings(flexsurv::flexsurvreg(survival::Surv(t, status) ~ 1, data=Data2, dist = "weibull", control = list(maxiter = 1000)))
  return(exp(c(survfit$coefficients[1],survfit$coefficients[2])))
}

#### Initial valures for the gamma parameters
initials.gamma <- function(Data){
  Data2 <- initials_data(Data)
  survfit <- suppressWarnings(flexsurv::flexsurvreg(survival::Surv(t, status) ~ 1, data=Data2, dist = "gamma"))
  return(c(exp(survfit$coefficients[1]),1/exp(survfit$coefficients[2])))
}

#### Initial valures for the lognormal parameters
initials.lnorm <- function(Data){
  Data2 <- initials_data(Data)
  survfit <- suppressWarnings(flexsurv::flexsurvreg(survival::Surv(t, status) ~ 1, data=Data2, dist = "lnorm"))
  return(c(survfit$coefficients[1],exp(survfit$coefficients[2])))
}

#### Initial valures for the log-logistic parameters
initials.llogis <- function(Data){
  Data2 <- initials_data(Data)
  survfit <- suppressWarnings(flexsurv::flexsurvreg(survival::Surv(t, status) ~ 1, data=Data2, dist = "llogis"))
  return(exp(c(survfit$coefficients[1],survfit$coefficients[2])))
}













