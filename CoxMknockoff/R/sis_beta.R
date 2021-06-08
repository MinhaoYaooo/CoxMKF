#' sis_beta
#'
#' @import ppcor
#' @description Sure independence screening based on M->Y pathway.
#'
#' @param X X is the vector of exposure variable.
#' @param M M is the matrix of mediators.
#' @param COV COV is the vector/matrix of covariate(s).
#' @param Y Y is the vector of outcome variable.
#'
#' @return a list containing variance, estimation and pval for beta

sis_beta <- function(X,M,Y,COV){
  s_beta <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(Y=Y, M = M[, j], X = X)
    } else {
      MX <- data.frame(Y=Y, M = M[, j], X = X, COV = COV)
    }
    fit <- coxph(Y ~., data = MX)
    s1 <- fit$var[1,1]                   #var for alpha
    s2 <- summary(fit)$coef[1]           #coefficients for alpha
    s3 <- summary(fit)$coef[1,5]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_beta)))
  beta_sis <- t(dat)
  colnames(beta_sis) = colnames(M)
  return(s_beta=beta_sis)
}
