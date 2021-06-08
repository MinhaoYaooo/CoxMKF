#' knockoff_surv
#'
#' @import corpcor
#' @import knockoff
#' @import ncvreg
#'
#' @description Knockoff method for survival outcome.
#'
#' @param X X is the vector of exposure variable.
#' @param Y Y is the vector of outcome variable.
#' @param M_SIS M_SIS is the matrix of mediators after sure independence screening.
#' @param COV COV is the vector/matrix of covariate(s).
#' @param penalty penalty is the type of penalty we use.
#' @param fdr fdr is the targeted false discovery rate that we want to control. A value between 0 and 1.
#' @param stat stat is the type of knockoff statistic we use.
#'
#' @return a vector of selected mediators.


knockoff_surv <- function(X, Y, M_SIS, COV, penalty=c("MCP", "SCAD", "lasso"), fdr, stat=c('beta', 'lambda')){

  mu = apply(M_SIS,2,mean)
  covariance = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(M_SIS,verbose=F)), nrow=ncol(M_SIS)))},
                        warning = function(w){}, error = function(e) {
                          stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                        }, finally = {})

  if(stat=='beta'){
    M_knockoff = knockoff::create.gaussian(M_SIS, mu, covariance)
  }else{
    M_knockoff = knockoff::create.fixed(M_SIS, method = 'equi')$Xk  # create the knockoff matrix
  }

  colnames(M_knockoff) <- paste0(colnames(M_SIS), '_K')

  XM <- cbind(M_SIS, M_knockoff, X)

  if (is.null(COV)) {
    fit <- ncvsurv(XM, Y,
                   penalty = penalty, nlambda = 5000,
                   penalty.factor = c(rep(1, 2*ncol(M_SIS)), 0))
  } else {
    COV <- data.frame(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
    conf.names <- colnames(COV)
    XM_COV <- cbind(XM, COV)
    fit <- ncvsurv(XM_COV, Y,
                   penalty = penalty, nlambda = 10000,
                   penalty.factor = c(rep(1, 2*ncol(M_SIS)), rep(0, 1 + ncol(COV))))
  }

  if(stat=='beta'){
    lam <- fit$lambda[which.min(BIC(fit))]
    # if(verbose) cat("lambda selected: ", lam, "\n")
    Coefficients <- coef(fit, lambda = lam)
    Z_stat <- Coefficients[1:(2*ncol(M_SIS))]
    W <- abs(Z_stat[1:ncol(M_SIS)]) - abs(Z_stat[(ncol(M_SIS)+1):(2*ncol(M_SIS))])
    names(W) <- colnames(M_SIS)[1:ncol(M_SIS)]
  }else{
    betas <- fit$beta
    Z_stat <- rep(0, 2*ncol(M_SIS))
    for (i in 1:(2*(ncol(M_SIS)))) {
      Z_stat[i] <- knockoff_stat(betas[i,], as.numeric(colnames(betas)))
    }
    W <- pmax(Z_stat[1:ncol(M_SIS)],
              Z_stat[(ncol(M_SIS)+1):(2*ncol(M_SIS))]) *
      sign(Z_stat[1:ncol(M_SIS)]-Z_stat[(ncol(M_SIS)+1):(2*ncol(M_SIS))])
    names(W) <- colnames(M_SIS)[1:ncol(M_SIS)]
  }


  thres <- threshold(W, fdr)

  ID_Knockoff <- which(W>=thres)

  return(ID_Knockoff)
}
