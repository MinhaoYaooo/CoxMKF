#' proportion
#'
#' @description Calculate the proportion of a given num, whicberh will be used in calculating the threshold.
#'
#' @param v v is a given number
#' @param W W is the vector of knockoff statistic.
#'
#' @return a value

proportion <- function(v, W){
  return(sum(W<=-v) / max(1, sum(W>=v)))
}
