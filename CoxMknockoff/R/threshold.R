#' threshold
#'
#' @description Get the threshold of knockoff statistic.
#'
#' @param W W is the vector of knockoff statistic.
#' @param fdr fdr is the targeted false discovery rate that we want to control. A value between 0 and 1.
#'
#' @return a value

threshold <- function(W, fdr){
  W_ <- unique(abs(W[W!=0]))
  W_ <- W_[order(W_)]
  props <- sapply(W_, proportion, W=W)
  return(W_[props<=fdr][1])
}
