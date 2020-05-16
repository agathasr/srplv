#' Selection criteria based on maximized log-likelihood function value
#'
#' \code{measures_loglik} calculates the AIC, AICc, BIC, HQIC, CAIC measures.
#' @param n  Number of systems.
#' @param npar  Number of the parameters.
#' @param loglik Maximized log-likelihood function value returned in \code{\link{EMalg}}.
#'
#' @return A list with the following components:
#'    \item{loglik}{The maximized log-likelihood function value.}
#'    \item{aic}{AIC (The Akaike Information Criterion).}
#'    \item{aicc}{AICc (The corrected Akaike Information Criterion).}
#'    \item{bic}{BIC (The Bayesian Information Criterion).}
#'    \item{hqic}{HQIC (The Hannan-Quinn Information Criterion).}
#'    \item{caic}{CAIC (The Consistent Akaike Information Criterion).}
#' @export
#'
measures_loglik <- function(n,npar,loglik){
  aic <- -2*loglik + 2*npar
  aicc <- aic + (2*npar*(npar+1))/(n-npar-1)
  bic <- -2*loglik + npar*log(n)
  hqic <- -2*loglik + 2*npar*log(log(n))
  caic <- -2*loglik + npar*(log(n)+1)
  return(list(loglik=loglik,aic=aic,aicc=aicc,bic=bic,hqic=hqic,caic=caic))
}


