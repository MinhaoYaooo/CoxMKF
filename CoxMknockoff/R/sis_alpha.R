#' sis_alpha
#'
#' @description Sure independence screening based on X->M pathway.
#'
#' @param X X is the vector of exposure variable.
#' @param M M is the matrix of mediators.
#' @param COV COV is the vector/matrix of covariate(s).
#'
#' @return a list containing variance, estimation and pval for alpha

sis_alpha <- function(X, M, COV){
  p <- ncol(M)
  s_alpha <- matrix(0, 3, p)
  for(j in 1:p){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s_alpha[1,j] <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s_alpha[2,j] <- summary(fit)$coef[2]           #coefficients for alpha
    s_alpha[3,j] <- summary(fit)$coef[2,4]         #p-value for alpha
  }
  colnames(s_alpha) = colnames(M)
  return(s_alpha = s_alpha)
}
