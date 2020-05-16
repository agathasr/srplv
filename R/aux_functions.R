#' GetTime
#' This function returns the current amount of spent time by the current R process (user time + system time). The unit is in minutes.
#'
.gettime <- function() {
  t <- proc.time()
  out <- (t[1]+t[2])/60
  names(out) <- NULL
  return(out)
}

#' This function calculates the mean of likelihood function evaluated at the vector par over L values.
lik.par <- function(par,Data,delta,L,distribution) {
  if(distribution=='lnorm'){
    par[2] <- exp(par[2])
  } else{
    par <- exp(par)
  }
  out <- 0
  for (l in 1:L) {
    out <- out+lik(par=par[1:2],Data=Data,delta=delta[,l],distribution)
  }
  out <- out/L
  return(out)
}


#' Function posteriori calculates the kernel of posterior density.

post.weibull <- function(par,Data,delta,distribution,prior1,prior2){
  par <- exp(par)
  prior.sh <- dgamma(par[1],shape=prior1[1],scale=prior1[2], log=TRUE)
  prior.scl <- dgamma(par[2],shape=prior2[1],scale=prior2[2], log=TRUE)
  post <- prior.sh + prior.scl + lik(par=par,Data=Data,delta=delta,distribution=distribution)
  return(post)
}

post.gamma <- function(par,Data,delta,distribution,prior1,prior2){
  par <- exp(par)
  prior.sh <- dgamma(par[1],shape=prior1[1],scale=prior1[2], log=TRUE)
  prior.scl <- dgamma(par[2],shape=prior2[1],scale=prior2[2], log=TRUE)
  post <- prior.sh + prior.scl + lik(par=par,Data=Data,delta=delta,distribution=distribution)
  return(post)
}

post.lnorm <- function(par,Data,delta,distribution,prior1,prior2){
  par[2] <- exp(par[2])
  prior.ml <- dnorm(par[1],mean=prior1[1],sd=sqrt(prior1[2]), log=TRUE)
  prior.vl <- dgamma(par[2],shape=prior2[1],scale=prior2[2], log=TRUE)
  post <- prior.ml + prior.vl + lik(par=par,Data=Data,delta=delta,distribution=distribution)
  return(post)
}

post.llogis <- function(par,Data,delta,distribution,prior1,prior2){
  par <- exp(par)
  prior.sh <- dgamma(par[1],shape=prior1[1],scale=prior1[2], log=TRUE)
  prior.scl <- dgamma(par[2],shape=prior2[1],scale=prior2[2], log=TRUE)
  post <- prior.sh + prior.scl + lik(par=par,Data=Data,delta=delta,distribution=distribution)
  return(post)
}


#### CPO function ##############

cpo.aux <- function(param,r,time,del,m,cens,distribution){
  if(distribution=='lnorm'){
    param[2] <- exp(param[2])
  } else{
    param <- exp(param)
  }
  out <- exp(eval(parse(text=paste('liksyst_', distribution, '(par=param,r=r,time=time,d=del,
                             m=m,cens=cens)', sep=''))))
  return(out)
}



