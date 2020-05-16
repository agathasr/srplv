#' Metropolis within Gibbs algorithm for srplv
#'
#' \code{MWG}  obtains MCMC sample from posterior distribution via MWG algorithm for  superposed renewal processes by means of latent variables.
#'
#' @param guess  The vector of initial values for the parameters. The default is (1,1).
#' @param Data  Data in the srplv format as \code{\link{data_to_srplv}}.
#' @param distribution The chose distribution: weibull, gamma, lnorm, llogis.
#' @param burn Burn-in sample.
#' @param jump Jump sample.
#' @param prior1 Hyperparameters for the prior distribution for parameter 1.
#' @param prior2 Hyperparameters for the prior distribution for parameter 2.
#' @param n.size Number of MCMC sample. The default is 1000.
#'
#' @return A list with the following components:
#'    \item{post.values}{The sample from the posterior distribution for each parameter (column) with size n.size.}
#'    \item{acceptance}{The acceptance rate.}
#'    \item{time}{The computational time.}
#'    \item{matrix.cpo}{Matrix with dimension n.sys by n.size containing Conditional Predictive Ordinate (CPO) values for n.sys systems and n.size posterior sample.}
#' @export
#'

MWG <- function(guess,Data,distribution,burn,jump,prior1,prior2,n.size=1000){
  if(!is.numeric(guess) | length(guess)!=2) {
    guess <- c(1,1)
    print("Not valid value for guess. Default values are used: (1,1).") }
  dist <- c("weibull","gamma","lnorm","llogis")
  if(distribution %in% dist==FALSE) {
    distribution='weibull'
    print("Distribution should be weibull, gamma, lnorm or llogis. The default is weibull.")
  }
  if(!is.list(Data)){
    print("Data sould be a list. See data_to_srplv for more detail.")
  }
  glue::glue('{format(Sys.time(), "%d/%b %H:%M:%S")} Lets run...')
  ini <- .gettime()
  if(distribution=='lnorm'){
    aux.par <- c(guess[1],log(guess[2]))
  } else{
    aux.par <- log(guess)
  }
  k <- 2
  while (k <= burn){
    aux.delta <- rdelta(par=aux.par,Data=Data,distribution=distribution)
    scl <- 0.4
    aux.metrop <- mcmc::metrop(eval(parse(text=paste('post.', distribution, sep=''))),
                               initial=aux.par,scale=scl,nbatch=1,Data=Data,delta=aux.delta,distribution=distribution,prior1=prior1,prior2=prior2)
    aux.par <- aux.metrop$batch
    k <- k+1
  }
  par.est <- matrix(NA,nrow=n.size,ncol=length(guess))
  k <- 2
  m <- 1
  ta <- 0
  ac <- 0
  matrix.cpo <- matrix(NA,nrow=Data$n,ncol=n.size)

  while (k <= n.size*jump){
    aux.delta <- rdelta(par=aux.par,Data=Data,distribution=distribution)
    scl <- ifelse(ac<0.3,0.1,
                  ifelse(ac>0.5,0.8,0.6))
    aux.metrop <- mcmc::metrop(eval(parse(text=paste('post.', distribution, sep=''))),
                               initial=aux.par,scale=scl,nbatch=1,Data=Data,delta=aux.delta,distribution=distribution,prior1=prior1,prior2=prior2)
    param.aux <- aux.metrop$batch
    if(k %% jump==0){
      par.est[m,] <- param.aux
      for(i in 1: Data$n){
        matrix.cpo[i,m] <- cpo.aux(param=param.aux,r=Data$r[i],time=Data$time[[i]],del=aux.delta[[i]],
                                   m=Data$n.comp,cens=Data$cens[i],distribution=distribution)
      }
      ta.aux <- ifelse(param.aux[1]==aux.par[1],0,1)
      ta <- ta+ta.aux
      ac <- ta/m
      m <- m+1
    }
    aux.par <- param.aux
    k <- k+1
  }
  spent <- .gettime() - ini
  print(glue::glue('Spent = {round(spent,3)} min.'))
  if(distribution=='lnorm'){
    out <- cbind(par.est[,1],exp(par.est[,2]))
  } else{
    out <- exp(par.est)
  }
  return(list(post.values=out,acceptance=ac,time=spent,matrix.cpo=matrix.cpo))
}

