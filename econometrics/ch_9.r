#
# Hayashi Ch. 9
#

# Monte Carlo Exercise

# Simulation 1

df_test <- function(yy){
  len <- length(yy)
  ols_mu <- lm(yy[-1]~yy[-len])
  t_mu <- as.double((coef(ols_mu)[2]-1) / summary(ols_mu)[[4]][2,2])
  t_trend <- c(1:(len-1))
  ols_tau <- lm(yy[-1]~yy[-len]+t_trend)
  t_tau <- as.double((coef(ols_tau)[2]-1) / summary(ols_tau)[[4]][2,2])
  return( c(t_mu,t_tau) )
}

sim1 <- function(nrep = 10^3, TT = 100, y_0 = 0, rho = 0.95){
  rr <- cumprod(rep(rho,TT))
  A <- diag(TT)
  for( j in 1:(TT-1))
    A[(j+1):TT, j] <- cumprod(rep(rho,TT-j))
  
  count_mu <- 0
  count_tau <- 0
  for( jj in 1:nrep){
    y <- rr * y_0 + A %*% rnorm(TT)
    ans <- df_test(c(y_0, y))
    if ( ans[1] < -2.86 ) count_mu <- count_mu+1
    if ( ans[2] < -3.41 ) count_tau <- count_tau+1
  }
  #cat(count_mu/nrep,count_tau/nrep,"\n")
return(c(count_mu/nrep,count_tau/nrep))
}

sim1()
#> sim1(10^4)
#[1] 0.1241 0.0953
#> sim1(3*10^4)
#[1] 0.1262333 0.0902000


# simulation 2

rm(list=ls())

adf_test <- function(y_t,nlag){
  len <- length(y_t)
  yy <- y_t[(nlag+2):len]
  nreg <- length(yy)
  rhs <- matrix(ncol=2+nlag,nrow=len-nlag-1)
  rhs[,1] <- 1
  rhs[,2] <- y_t[(nlag+1):(len-1)]
  dif_y <- diff(y_t[1:(len-1)])
  for (jj in 1:nlag)
    rhs[,jj+2] <- dif_y[(nlag+1-jj):(len-1-jj)]
  ols <- lm(yy~rhs-1)
  coefs <- coef(ols)
  rho <- as.double(nreg*(coefs[2]-1)/(1-sum(coefs[3:(nlag+2)])))
  tau <- as.double((coefs[2]-1) / summary(ols)[[4]][2,2])
  return( c(rho,tau) )
}

sim2 <- function(nrep=10^3, TT = 100, y_0 = 0, e_0 = 0, nlag=4, theta = -0.8){
  B <- diag(TT)
  C <- diag(TT)
  for( j in 1:TT){
    B[j:TT, j] <- 1
    C[j:TT, j] <- theta
  }
  
  count_rho <- 0
  count_t <- 0
  for( jj in 1:nrep){
    e_t <- rnorm(TT)
    y <- B %*% e_t + C %*% c(e_0,e_t[-TT])
    ans <- adf_test(c(y_0, y),nlag)
    if ( ans[1] < -14.1 ) count_rho <- count_rho+1
    if ( ans[2] < -2.86 ) count_t <- count_t+1
  }
  #cat(count_rho/nrep,count_t/nrep,"\n")
return(c(count_rho/nrep,count_t/nrep))
}


#> sim1(10^4)
#[1] 0.1241 0.0953


sim2()
#> sim2(10^5)
#[1] 0.49702 0.28952

# Empirical Exercise

library(urca)
library(tseries)

data_lt <- read.table("lt.asc", header=F)

S_t <- data_lt[[2]]
P_t <- data_lt[[3]]
P_t_st <- data_lt[[4]]
nobs <- length(S_t)


z_t <- log(S_t) - log(P_t) + log(P_t_st)

z_t_ts <- ts(z_t, start=1791)

# (a)

pmax <- 14
zz <- z_t[(pmax+2):nobs]
nreg <- length(zz)
zmat <- matrix(1, nrow=nobs-pmax-1, ncol=pmax+2)
zmat[,2] <- z_t[(pmax+1):(nobs-1)]
diffz <- diff(z_t, 1)

for (j in 1:pmax)
  zmat[,j+2] <- diffz[(pmax-j+1):(nobs-1-j)]

lmlist <- vector("list", pmax+1)

AIC_z <- c()
BIC_z <- c()
for( j in 0:pmax){
  lmlist[[j+1]] <- lm(zz ~ zmat[,1:(j+2)]-1)
  SSR <- deviance(lmlist[[j+1]])
  AIC_z[j+1] <- log(SSR/nreg) + (j+2)*2/nreg
  BIC_z[j+1] <- log(SSR/nreg) + (j+2)*log(nreg)/nreg
  }

summary(lmlist[[15]])
coefs <- coef(lmlist[[15]])
#adf_rho <- as.double(nreg*(coefs[2]-1) /( 1 - sum(coefs[3:16]) ))
adf_tau <- as.double((coefs[2]-1) / summary(lmlist[[15]])[[4]][2,2])

which.min(AIC_z)
which.min(BIC_z)

summary(ur.df(z_t, type = c("drift"), lags = 14))
summary(ur.df(z_t, type = c("drift"), lags = 0))
adf.test( z_t, k=14 )

# (b)

z_b <- window(z_t_ts, 1974)
nreg_b <- length(z_b)-1
ols_b <- lm(z_b[-1] ~ z_b[-(nreg_b+1)] )

adf_tau_b <- as.double((coef(ols_b)[2]-1) / summary(ols_b)[[4]][2,2])
adf_tau_b

# (c)

z_c <- window(z_t_ts, 1870, 1913)
nreg_c <- length(z_c)-1
ols_c <- lm(z_c[-1] ~ z_c[-(nreg_c+1)] )

adf_tau_c <- as.double((coef(ols_c)[2]-1) / summary(ols_c)[[4]][2,2])
adf_tau_c

# (d)

KK <- 2
SSR_r <- deviance(lm(z_t[-1] ~ z_t[-nobs] ))

z_1 <- window(z_t_ts, 1791, 1974)
nreg_1 <- length(z_1)-1
SSR_1 <- deviance(lm(z_1[-1] ~ z_1[-(nreg_1+1)] ))

z_2 <- window(z_t_ts, 1974)
nreg_2 <- length(z_2)-1
SSR_2 <- deviance(lm(z_2[-1] ~ z_2[-(nreg_2+1)] ))

F_d <- (SSR_r - SSR_1 - SSR_2)/KK/((SSR_1+SSR_2)/(nobs-1-2*KK))
F_d * KK
1 - pchisq(F_d*KK,KK)

# (e)

TT <- length(z_c)
rho_bar <- 1 - 7/TT
lhs <- c()
lhs[1] <- z_c[1]
lhs[2:TT] <- z_c[2:TT] - rho_bar * z_c[1:(TT-1)]
rhs <- c()
rhs[1] <- 1
rhs[2:TT] <- 1 - rho_bar

gls <- lm(lhs~rhs-1)
b_bls <- coef(gls)

y_til <- z_c - b_bls

pmax_e <- 9
zz_e <- y_til[(pmax_e+2):TT]
zz_lag <- y_til[(pmax_e+1):(TT-1)]
nreg_e <- length(zz_e)

diffz_e <- diff(y_til, 1)

zmat_e <- matrix(ncol=pmax_e +1, nrow=nreg_e)
zmat_e[,1] <- zz_lag
for (j in 1:pmax_e)
  zmat_e[,j+1] <- diffz_e[(pmax_e-j+1):(TT-1-j)]

BIC_e <- c()
for (j in 0:pmax_e)
  {
  SSR <- deviance(lm(zz_e ~ zmat_e[,1:(j+1)]-1))
  BIC_e[j+1] <- log(SSR/nreg_e) + (j+2)*log(nreg_e)/nreg_e
  }
which.min(BIC_e)

summary(ur.df(y_til, type = c("none"), lags = 0))





