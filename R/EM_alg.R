#' EM algorithm for srplv
#'
#' \code{EMalg} returns the maximum likelihood estimates obtained by EM algorithm.
#'
#' @param guess  The initial values for the parameters. The default is (1,1).
#' @param Data  Data in the srplv format as \code{\link{data_to_srplv}}.
#' @param distribution The chose distribution: weibull, gamma, lnorm, llogis.
#' @param L  Number of Monte Carlo simulation.
#'
#' @return A list with the following components:
#'    \item{values}{The values of parameters in all EM algorithm iterations.}
#'    \item{est.par}{The maximum likelihood estimates of the parameters.}
#'    \item{max.logLik}{The maximum likelihood value.}
#'    \item{n.it}{The number of EM iterations.}
#'    \item{time}{The computational time.}
#'
#' @export
#'
EMalg <- function(guess,Data,distribution,L){
  if(!is.numeric(guess) | length(guess)!=2) {
    guess <- c(1,1)
    print("Not valid value for guess. Default values are used: (0,0).") }
  dist <- c("weibull","gamma","lnorm","llogis")
  if(distribution %in% dist==FALSE) {
    distribution='weibull'
    print("Distribution should be weibull, gamma, lnorm or llogis. The default is weibull.")
  }
  if(L<=0){
    L <- 1000
    print("Not valid value of L. The default is 1000.")
  }
  if(!is.list(Data)){
    print("Data sould be a list. See data_to_srplv documentation for more detail.")
  }
  if(distribution=='lnorm'){
   guess1 <- c(guess[1],log(guess[2]))
  } else{
   guess1 <- log(guess)
  }
  est.param <- matrix(0,ncol=2,nrow=200)
  est.param[1,] <- guess1
  glue::glue('{format(Sys.time(), "%d/%b %H:%M:%S")} Lets run...')
  ini <- .gettime()
  it <- 1
  aux.delta <- replicate(L,rdelta(par=guess,Data=Data,distribution),simplify="array")
  res <- optim(guess,lik.par,Data=Data,delta=aux.delta,L=L,distribution=distribution,control=list(fnscale=-1))
  est.param[it+1,] <- res$par
  spent <- .gettime() - ini
  print(glue::glue('Iteration {it}: Distribution={distribution}; Parameters = {round(exp(res$par[1]),3)};{round(exp(res$par[2]),3)}; Spent = {round(spent,3)} min.'))
  it <- 2
  aux.delta <- replicate(L,rdelta(par=res$par,Data=Data,distribution),simplify="array")
  aux <- optim(res$par,lik.par,Data=Data,delta=aux.delta,L=L,distribution=distribution,control=list(fnscale=-1))
  est.param[it+1,] <- aux$par
  spent <- .gettime() - ini
  print(glue::glue('Iteration {it}: Distribution={distribution}; Parameters = {round(exp(aux$par[1]),3)};{round(exp(aux$par[2]),3)}; Spent = {round(spent,3)} min.'))
  while (((abs(aux$par[1] - res$par[1]) > 0.0001) | (abs(aux$par[2] - res$par[2]) > 0.0001)) & (it<199)) {
    res <- aux
    it <- it+1
    aux.delta <- replicate(L,rdelta(par=res$par,Data=Data,distribution),simplify="array")
    aux <- optim(res$par,lik.par,Data=Data,delta=aux.delta,L=L,distribution=distribution,control=list(fnscale=-1))
    est.param[it+1,] <- aux$par
    spent <- .gettime() - ini
    print(glue::glue('Iteration {it}: Distribution={distribution}; Parameters = {round(exp(aux$par[1]),3)};{round(exp(aux$par[2]),3)}; Spent = {round(spent,3)} min.'))
  }
  if(distribution=='lnorm'){
    out.est <- cbind(est.param[1:(it+1),1],exp(est.param[1:(it+1),2]))
    out.par <- c(aux$par[1],exp(aux$par[2]))
  } else{
    out.est <- exp(est.param[1:(it+1),])
    out.par <- exp(aux$par)
  }
  delta.ger <- replicate(L,rdelta(par=aux$par,Data=Data,distribution),simplify="array")
  max.logLik <- lik.par(par=c(aux$par[1],aux$par[2]),Data=Data,delta=delta.ger,L=L,distribution=distribution)
  return(list(values=out.est,est.par=out.par,max.logLik=max.logLik,n.it=it,time=spent))
}



