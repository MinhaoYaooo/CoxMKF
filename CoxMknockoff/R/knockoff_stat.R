#' knockoff_stat
#'
#' @description Calculate the knockoff statistic based on the result of penalized regression.
#'
#' @param lambda lambda is the vector of penalties.
#' @param beta beta is the vector of coefficients corresponding to lambda.
#'
#' @return a vector

knockoff_stat <- function(beta, lambda){
  n <- length(beta)
  Z_stat <- lambda[1]
  for (i in 2:n) {
    if((beta[i-1]==0)&(beta[i]!=0)){
      Z_stat <- lambda[i]
    }
  }
  if(beta[n]==0){
    Z_stat <- lambda[n]
  }
  return(Z_stat)
}
