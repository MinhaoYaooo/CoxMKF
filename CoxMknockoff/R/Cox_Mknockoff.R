#' Cox_Mknockoff
#'
#' @import corpcor
#' @import knockoff
#' @import ncvreg
#' @import ppcor
#' @import MultiMed
#'
#' @description Apply knockoff methods into high-dimensional mediation analysis.
#'
#' @param X X is the vector of exposure variable.
#' @param Y Y is the vector of outcome variable.
#' @param M M is the matrix of mediators.
#' @param COV COV is the vector/matrix of covariate(s).
#' @param k k is an integer to adjust the number of mediators in SIS step.
#' @param penalty penalty is the type of penalty we use.
#' @param fdr fdr is the targeted false discovery rate that we want to control. A value between 0 and 1.
#' @param stat stat is the type of knockoff statistic we use.
#' @param knockoff knockoff is a boolean value to control whether we use knockoff method or not.
#' @param topN topN is an integer. If given, then we directly choose topN mediators in SIS step based on the statistics.
#' @param verbose if TRUE, then details will be printed
#'
#' @return a list of outcomes containing the selected mediators and the corresponding results.
#'
#'

Cox_Mknockoff <- function(X, Y, M, COV,k, fdr,
                          penalty = c("MCP", "SCAD", "lasso"),
                          stat = c('beta', 'lambda'),
                          knockoff = c(TRUE, FALSE),
                          topN = NULL,
                          verbose = TRUE) {

  penalty <- match.arg(penalty)

  n <- nrow(M)
  p <- ncol(M)

  if (is.null(topN)) {
    d <- ceiling(k*n/log(n))  #the top d mediators that associated with outcome
  } else {
    d <- topN
  }

  if(verbose) message("Step 1: Prelimenary Screening...", "     (", Sys.time(), ")")

  # SIS for alpha
  alpha_s <- sis_alpha(X, M, COV)
  SIS_alpha <- alpha_s[3,]
  SIS_alpha_sort <- sort(SIS_alpha)

  # SIS for beta
  beta_s <- sis_beta(X=X, M=M, Y=Y, COV=COV)

  gamma_s <- alpha_s[2,] * beta_s[2, ]
  gamma_abs_sort <- sort(abs(gamma_s), decreasing = T)

  ID_SIS <- which(abs(gamma_s)>=gamma_abs_sort[d])

  M_SIS <- M[, ID_SIS]

  if(verbose) cat("Top", length(ID_SIS), "mediators selected (ranked from most to least significant): ", colnames(M_SIS), "\n")


  if(verbose) message("Step 2: Penalized Variable Selection (", penalty, ") ...", "   s  (",
                      Sys.time(), ")")

  if(knockoff==TRUE){
    ID_Knockoff <- knockoff_surv(X=X, Y=Y, M_SIS=M_SIS, COV=COV, penalty = penalty, fdr = fdr, stat = stat)
  }else{
    XM <- cbind(M_SIS, X)

    if (is.null(COV)) {
      fit <- ncvsurv(XM, Y,
                     penalty = penalty, nlambda = 5000,
                     penalty.factor = c(rep(1, ncol(M_SIS)), 0))
    } else {
      COV <- data.frame(COV)
      COV <- data.frame(model.matrix(~., COV))[, -1,drop=F]
      conf.names <- colnames(COV)
      XM_COV <- cbind(XM, COV)
      fit <- ncvsurv(XM_COV, Y,
                     penalty = penalty, nlambda = 10000,
                     penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(COV))))
    }

    lam <- fit$lambda[which.min(BIC(fit))]
    Coefficients <- coef(fit, lambda = lam)
    est <- Coefficients[1:length(ID_SIS)]
    ID_Knockoff <- which(est != 0)
  }


  ID_Knockoff <- ID_SIS[ID_Knockoff]

  if(verbose) cat("Non-zero", penalty, "beta estimate(s) of mediator(s) found: ", colnames(M)[ID_Knockoff], "\n")

  XM <- cbind(M[,ID_Knockoff], X)
  C_M <- colnames(XM)

  Knockoff_M <- colnames(M)[ID_Knockoff]


  if(verbose) message("Step 3: The adjusted Sobel significance test ...",
                      "     (", Sys.time(), ")")

  alpha_s <- as.matrix(alpha_s[, Knockoff_M])

  alpha_est <- alpha_s[2, ]   #  the estimator for alpha
  var_alpha <- alpha_s[1, ]

  # true alpha and beta
  #beta_t <- beta[ID_Knockoff]
  #alpha_t <- alpha[ID_Knockoff]
  #ab_true <- alpha_t * beta_t

  if (is.null(COV)) {
    YMX <- data.frame(Y = Y, M[, ID_Knockoff, drop = FALSE], X = X)
  } else {
    YMX <- data.frame(Y = Y, M[, ID_Knockoff, drop = FALSE], X = X, COV = COV)
  }

  cox_model <- coxph(Y ~ ., data = YMX)
  beta_est <- summary(cox_model)$coefficients[1: length(ID_Knockoff)]     #the estimator of beta
  DE <- summary(cox_model)$coefficients[(length(ID_Knockoff)+1), 1]
  DE <- exp(DE)
  if(length(ID_Knockoff)==1){
    var_beta <- cox_model$var[1:length(ID_Knockoff),1:length(ID_Knockoff)]
  }else{
    var_beta <- diag(cox_model$var[1:length(ID_Knockoff),1:length(ID_Knockoff)])
  }
  pMY <- summary(cox_model)$coefficients[1:length(ID_Knockoff),5]

  pEM <- SIS_alpha[ID_Knockoff]

  multi_test <- medTest.SBMH(pEM = pEM, pMY = pMY, MCP.type = "FDR")

  ab_est <- alpha_est * beta_est   # the estimator of alpha*beta

  # var(alpha*beta)
  var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2)

  # confidence interval
  conf_low <- ab_est - 1.96 * sqrt(var_ab)
  conf_up <- ab_est + 1.96 * sqrt(var_ab)

  results <- data.frame(alpha = alpha_est, beta = beta_est,
                        `alpha_est*beta_est` = ab_est, #`alpha_t*beta_t` = ab_true,
                        conf_low=conf_low, conf_up=conf_up,
                        P_multimed=multi_test, var_ab=var_ab,
                        var_alpha=var_alpha, var_beta=var_beta,
                        DE, check.names = FALSE)

  if(verbose) message("Done!", "     (", Sys.time(), ")")

  return(list(SIS=C_M, Selected=Knockoff_M, Result=results))
}




