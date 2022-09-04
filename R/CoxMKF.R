#' @title Screening Based on X->M with FDR control
#' @description The first step of CoxMKF, which conducts preliminary screening with the p-values from X to M with FDR control.
#' @param X a vector of the exposure.
#' @param M a matrix of the mediators.
#' @param COV optional, a matrix of the potential covariates.
#' @param intercept binary, whether an intercept should be included in the regression M~X, the default value is TRUE.
#' @param fdr the pre-defined FDR level in the screening step, the default value is 0.2.
#'
#' @return A list of the screening results including the IDs, pvalues and coefficients.
#' @export

screen_FDR <- function(X, M, COV, intercept=TRUE, fdr=0.2){
  p <- ncol(M)
  p_alpha <- rep(NA, p)
  coef_alpha <- rep(NA, p)
  for(j in 1:p){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    if(intercept){
      fit <- stats::lm(M ~1+., data = MX)
      p_alpha[j] <- summary(fit)$coef[2,4]   #p-value for alpha
      coef_alpha[j] <- summary(fit)$coef[2,1]
    }else{
      fit <- stats::lm(M ~0+., data = MX)
      p_alpha[j] <- summary(fit)$coef[1,4]   #p-value for alpha
      coef_alpha[j] <- summary(fit)$coef[1,1]
    }
  }
  names(p_alpha) <- colnames(M); names(coef_alpha) <- colnames(M)
  p_alpha.adjust <- p.adjust(p_alpha, method = 'fdr')
  ID_screen <- which(p_alpha.adjust<=fdr)
  return(list(ID_screen=ID_screen, pval=p_alpha, coef=coef_alpha))
}

#' @title Knockoff for Survival Response
#' @describeIn Generate Knockoff copies for survival outcome and compute the corresponding statistics.
#' @param X a vector of the exposure.
#' @param Y a vector of the outcome.
#' @param MS a matrix of the mediators after the pre-liminary screening.
#' @param COV optional, a matrix of the potential covariates.
#' @param penalty the penalty used in the cox regression.
#' @param q the pre-defined FDR level in the knockoff procedure.
#' @param gamma the pre-defined quantile level.
#' @param n_bootstraps the pre-defined number of knockoffs in the AKO.
#'
#' @return A vector of mediator IDs selected by AKO.
#' @export

knockoff_surv <- function(X, Y, MS, COV, penalty=c('MCP','lasso','SCAD'),  q=0.1, gamma = 0.05, n_bootstraps = 25){

  pS <- ncol(MS)
  pvals = matrix(0, pS, n_bootstraps)

  for (i in 1:n_bootstraps) {
    #print(i)

    mu <- apply(MS,2,mean)
    Sig <- matrix(as.numeric(corpcor::cov.shrink(MS,verbose=F)), nrow=ncol(MS))

    MK = knockoff::create.gaussian(MS, mu, Sig)
    colnames(MK) <- paste0(colnames(MS))

    XM <- cbind(MS, MK, X)

    if (is.null(COV)) {
      fit <- ncvreg::ncvsurv(XM, Y,
                     penalty = penalty, nlambda = 5000,
                     penalty.factor = c(rep(1, 2*pS), 0))
    }else {
      XM_COV <- cbind(XM, COV)
      fit <-  ncvreg::ncvsurv(XM_COV, Y,
                     penalty = penalty, nlambda = 1000,
                     penalty.factor = c(rep(1, 2*pS), rep(0, 1 + ncol(COV))))
    }

    lam <- fit$lambda[which.min(BIC(fit))]
    Coefficients <- coef(fit, lambda = lam)

    w = c()

    for (m in 1:pS) {
      w = c(w, abs(Coefficients[m])-abs(Coefficients[m + pS]))
    }

    pvals[,i] = empirical_pval(w, offset = 0)

  }
  aggregated_pval = data.frame()
  aggregated_pval[1, 1:pS] = apply(pvals, 1, quantile_aggregation, gamma=gamma, gamma_min=gamma_min, adaptive=FALSE)
  colnames(aggregated_pval) = paste0(colnames(MK))

  threshold = fdr_threshold(aggregated_pval[1,], fdr=q, method='bhq',
                            reshaping_function=NULL)

  ID_Knockoff = c()

  for(m in 1:pS){if(aggregated_pval[1,m]<=threshold){ID_Knockoff = c(ID_Knockoff, paste0(colnames(aggregated_pval)[m]))}}

  return(ID_Knockoff)
}


#' @title High Dimensional Mediation Analysis using Multiple Knockoffs
#' @description The main function of CoxMKF which combines the screening step, knockoff step and estimation step.
#' @param X a vector of the exposure.
#' @param Y a vector of the outcome.
#' @param M a matrix of the mediators.
#' @param COV optional, a matrix of the potential covariates.
#' @param penalty the penalty used in the cox regression.
#' @param intercept binary, whether an intercept should be included in the regression M~X, the default value is TRUE.
#' @param q1 the pre-defined FDR level in the screening step.
#' @param q2 the pre-defined FDR level in the knockoff procedure.
#' @param gamma the pre-defined quantile level.
#' @param n_bootstraps the pre-defined number of knockoffs in the AKO.
#'
#' @return A list of the CoxMKF results including the IDs of the screening step, the IDs of the selection step and the estimation for the effects.
#' @export


CoxMKF <- function(X, Y, M, COV, penalty=c('MCP', 'lasso', 'SCAD'), intercept=TRUE, q1=0.2,
                   q2=0.1, gamma = 0.05, n_bootstraps = 25){

  n <- nrow(M); p <- ncol(M)

  ##### STEP 1: Preliminary Screening #####

  screen.results <- screen_FDR(X, M, COV, intercept = intercept, fdr=q1)

  MS <- M[, screen.results$ID_screen]

  ##### STEP 2: Selection with Knockoff #####

  ID_Knockoff <- knockoff_surv(X, Y, MS, COV, penalty = penalty,q=q2, gamma = gamma, n_bootstraps = n_bootstraps)

  ID_Knockoff <- screen.results$ID_screen[ID_Knockoff]

  ##### STEP 3: Estimation of Effects #####

  alpha_est <- screen.results$coef[ID_Knockoff]

  DATA <- data.frame(Y = Y, M[, ID_Knockoff], X = X, COV = COV)
  cox_model <- survival::coxph(Y ~ ., data = DATA)

  beta_est <- summary(cox_model)$coefficients[1: length(ID_Knockoff)]

  result <- data.frame(alpha=alpha_est, beta=beta_est, `alpha*beta` = alpha_est*beta_est)

  return(list(first.step=names(screen.results$ID_screen),
              second.step=names(ID_Knockoff),
              estimates=result))
}




########## Useful Functions for AKO ##########

bhy_threshold = function(pvals,
                         reshaping_function = NULL, fdr = 0.1)
{
  # """Benjamini-Hochberg-Yekutieli procedure for controlling FDR, with input
  #   shape function. Reference: Ramdas et al (2017)
  # """
  n_features = length(pvals)
  p_vals_sorted = sort(pvals)
  selected_index = 2 * n_features

  # Default value for reshaping function -- defined in
  # Benjamini & Yekutieli (2001)

  if (is.null(reshaping_function)==TRUE)
  {
    temp = seq(n_features,1)
    sum_inverse = 0
    for (i in temp) {
      sum_inverse = sum_inverse + 1 / i
    }

    return (bhq_threshold(pvals, fdr / sum_inverse))
  }
  else{
    for (i in seq(n_features - 1, 0, -1)){
      if (p_vals_sorted[i] <= fdr * reshaping_function(i) / n_features){
        selected_index = i
        break
      }

    }

    if (selected_index <= n_features){
      return (p_vals_sorted[selected_index])
    }
    else{
      return ('-1.0')
    }

  }

}



bhq_threshold = function(pvals, fdr=0.1){
  #   """Standard Benjamini-Hochberg for controlling False discovery rate
  # """
  n_features = length(pvals)
  pvals_sorted = sort(pvals)
  selected_index = 2 * n_features
  for (i in seq(n_features, 1, -1)){
    if (pvals_sorted[i] <= (fdr * i / n_features)){
      selected_index = i
      break
    }
  }

  if (selected_index <= n_features){
    return (pvals_sorted[selected_index])
  }

  else{
    return ('-1.0')
  }
}


fdr_threshold = function(pvals, fdr=0.1, method='bhq', reshaping_function=NULL){
  if (method == 'bhq'){
    # pvals_bhq = as.vector(unlist(pvals))
    return (bhq_threshold(pvals, fdr=fdr))
  }
  else{
    if(method == 'bhy'){
      # pvals_bhy = as.vector(unlist(pvals))
      return( bhy_threshold(
        pvals, fdr=fdr, reshaping_function=reshaping_function))
    }
    else{
      return('{} is not support FDR control method')
    }
  }
}


empirical_pval = function(test_score, offset = 1){
  pvals = c()
  n_features = length(test_score)
  if (offset !=0 && offset!=1){
    return("'offset' must be either 0 or 1")
  }
  else{
    test_score_inv = -test_score
    for (i in 1:n_features){
      if (test_score[i] <= 0){
        pvals = c(pvals, 1)
      }
      else{
        # pvals = c(pvals,(offset+sum(test_score_inv >= test_score[i]))/n_features)
        pvals = c(pvals,(offset+sum(test_score_inv[i] >= test_score))/n_features)
      }
    }
  }
  return (pvals)
}


fixed_quantile_aggregation = function(pvals, gamma = 0.05){
  #   """Quantile aggregation function based on Meinshausen et al (2008)
  # Parameters
  # ----------
  # pvals : 2D ndarray (n_bootstrap, n_test)
  # p-value (adjusted)
  # gamma : float
  # Percentile value used for aggregation.
  # Returns
  # -------
  # 1D ndarray (n_tests, )
  # Vector of aggregated p-value
  # """
  # pvals_fixed = as.vector(unlist(pvals))
  converted_score = (1 / gamma) *  quantile(pvals, gamma)

  return (min(1, converted_score))
}




adaptive_quantile_aggregation = function(pvals, gamma_min=0.05){
  #   """adaptive version of the quantile aggregation method, Meinshausen et al.
  # (2008)"""
  gammas = seq(gamma_min, 1.05, 0.05)
  list_Q = matrix(0,nrow = length(gammas), ncol = length(gammas))
  for (gamma in gammas) {
    list_Q = c(list_Q, fixed_quantile_aggregation(pvals, gamma))

  }

  return (min(1, (1 - log(gamma_min)) * min(list_Q)))

}




quantile_aggregation = function(pvals, gamma=0.5, gamma_min=0.00001, adaptive=FALSE){
  if (adaptive == TRUE){
    return (adaptive_quantile_aggregation(pvals, gamma_min))
  }

  else{
    return (fixed_quantile_aggregation(pvals, gamma))
  }

}






mfdr_tpp_knockoff = function(q, w, t, beta_true, beta_est){
  Sa = c()
  for(m in 1:length(w)){if(w[m]<=t){Sa = c(Sa, m)}}

  true_beta = which(beta_true!=0)
  Ec = which(beta_true == 0)
  mFDR_k = length(intersect(Sa, Ec))/(length(Sa)+1/q)

  #compute TPP with knockoffs

  tpp_k= length(intersect(Sa, true_beta))
  TPP_k = tpp_k/max(1,length(true_beta))

  return(list(mFDR_withknockoff = mFDR_k, TPP_withknockoff = TPP_k, l.sa = length(Sa)))
}


#compute baseline mFDR and TPP
mfdr_tpp_base = function(q, beta_true, beta_est){
  true_beta = which(beta_true!=0)
  Ec = which(beta_true ==0)

  Sa_b = which(beta_est!=0)
  mFDR_b = length(intersect(Sa_b, Ec))/(length(Sa_b)+1/q)

  #compute TPP
  t_base = length(intersect(Sa_b, true_beta))
  TPP_base = t_base/max(1,length(true_beta))

  return(list(mFDR_baseline = mFDR_b, TPP_baseline = TPP_base))
}



