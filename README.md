# CoxMKF
`CoxMKF` is an R package to apply the aggregated knockoffs to high dimensional mediation analysis with a survival response.

# Examples
First, we generate the data:
```
n <- 500                                   #number of samples
p <- 1000                                  #number of mediators
alpha=rep(0,p)                             #coefficients (mediator~exposure)
beta=rep(0,p)                              #coefficients (outcome~mediators)
alpha[1:12] <- c(0.55,0.45,-0.4,-0.45,0.5,0.6,-0.4,-0.46,-0.4,0.5,0,0)
beta[1:12] <- c(0.52,0.45,0.4,0.4,-0.54,-0.6,-0.4,-0.5,0,0,0.4,-0.8)
X <- t(t(rbinom(n, 1, 0.6)))               #exposure
Z1 <- t(t(rbinom(n, 1, 0.3)))              #covariates Z1
theta1 <- 0.3                              #coefficients(Z1-->M)
Z2 <- t(t(runif(n, 0, 1)))                 #covariates Z2
theta2 <- 0.2                              #coefficients(Z2-->M)
Z <- cbind(Z1, Z2)
phi <- c(0.3, -0.2)                        #coefficients(covariates-->outcome)
ck <- t(runif(p, 0, 1))
M <- matrix(0, n, p)                       #mediators
for(i in 1:n){
  e <- rnorm(p, sd = 1)
  M[i,] <- ck+X[i]*alpha+Z[i,1]*theta1+Z[i,2]*theta2+e
}
colnames(M) <- paste0("M", 1:ncol(M))
haz <- 0.5*exp(0.5*X+0.3*Z[,1]-0.2*Z[,2]+M%*%beta)   #baseline hazard function lambda0 <- 0.5
ft <- rexp(n, haz)
ct <- rexp(n, 0.7)                         #censoring time
time <- pmin(ft, ct)                       #observed time
status <- as.numeric(ft <= ct)             #censoring indicator
Y <- survival::Surv(time, status)
COV <- Z
```
Then we apply `CoxMKF` to select the mediators:
```
results <- CoxMKF(X, Y, M, COV, penalty = 'MCP', q2 = q, gamma = gamma, n_bootstraps = n_bootstraps)
```

