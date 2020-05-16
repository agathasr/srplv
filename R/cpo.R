#' CPO - Conditional Predictive Ordinate - for srplv
#'
#' \code{cpo_fc} calculates the sum of logarithm of CPO.
#'
#' @param matrix.cpo  Matrix with dimension n.sys by n.size containing CPO values for n.sys systems and n.size posterior sample. Result from \code{\link{MWG}}.
#'
#' @return sum of log of CPO.
#' @export
#'
cpo_fc <- function(matrix.cpo){
  cpo.i <- rep(NA,dim(matrix.cpo)[1])

  for (i in 1:dim(matrix.cpo)[1]){
    out <- (1/dim(matrix.cpo)[2])*sum(1/matrix.cpo[i,])
    cpo.i[i] <- 1/out
  }
  CPO.comp <- sum(log(cpo.i))
  return(CPO.comp)
}
