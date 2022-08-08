#
# Hayashi Ch. 10
#

# Empirical Exercise

mpyr <- read.table("mpyr.asc",header=F)

logm1 <- mpyr[[1]]
logp <- mpyr[[2]]
logy <- mpyr[[3]]
rr <- mpyr[[4]]
nobs <- length(rr)


# (a)

library(urca)
summary(ur.df(y=rr, lags=2, type='drift'))
summary(ur.df(y=rr, lags=4, type='drift'))

summary(ur.df(y=logm1-logp, lags=2, type='trend'))
summary(ur.df(y=logm1-logp, lags=4, type='trend'))

summary(ur.df(y=logy, lags=2, type='trend'))
summary(ur.df(y=logy, lags=4, type='trend'))


# (b)

ols_a1 <- lm(logm1-logp ~ logy + rr  )
summary(ols_a1)
x <- resid(ols_a1)

dx_l <- diff(x)[-1]
x_r <- x[-c(1,length(x))]
d_x1 <- diff(x)[-(length(x)-1)]
ols_a2 <- lm(dx_l ~ x_r + d_x1 - 1)
summary(ols_a2)
t_value <- summary(ols_a2)[[4]][1,3]
print(t_value)

# (c)

sols <- lm(logm1[4:88]-logp[4:88] ~ logy[4:88] + rr[4:88]  )
summary(sols)

dy <- diff(logy)
dr <- diff(rr)

rhs <- matrix(1,ncol=13,nrow=85)
rhs[,2] <- logy[4:88]
rhs[,3] <- rr[4:88]
rhs[,4] <- dy[1:85]
rhs[,5] <- dy[2:86]
rhs[,6] <- dy[3:87]
rhs[,7] <- dy[4:88]
rhs[,8] <- dy[5:89]
rhs[,9] <- dr[1:85]
rhs[,10] <- dr[2:86]
rhs[,11] <- dr[3:87]
rhs[,12] <- dr[4:88]
rhs[,13] <- dr[5:89]
dols <- lm(logm1[4:88]-logp[4:88] ~ rhs -1)
summary(dols)

# (d)

D_t <- c(rep(0,43),rep(1,42))

XX <- matrix(1,ncol=16,nrow=85)
XX[,2:3]<-rhs[,2:3]
XX[,4] <- D_t
XX[,5] <- logy[4:88]*D_t
XX[,6] <- rr[4:88]*D_t
XX[,7:16]<-rhs[,4:13]

dols_d <- lm(logm1[4:88]-logp[4:88] ~ XX -1)
summary(dols_d)
std <- sqrt(deviance(dols_d)/(85-16))
coef_dols <- coef(dols_d)
v_hat <- resid(dols_d)

ols_d <- lm(v_hat[-c(1,2)]~v_hat[c(-1,-85)]+v_hat[c(-84,-85)]-1)
summary(ols_d)
coefs <- coef(ols_d)
sig <- sqrt(deviance(ols_d)/(83-16-2))
lam_hat <- sig/(1-sum(coefs))

coef_dols[2:6]
sqrt(diag(solve(t(XX)%*%XX))[2:6])*lam_hat


# Wald test (–¢Š®¬)

nobs <- 85
RR <-  matrix(c(
c(0,0,0,1,0,0,rep(0,10)),
c(0,0,0,0,1,0,rep(0,10)),
c(0,0,0,0,0,1,rep(0,10))),nrow=3,byrow=T)
a_b <- RR %*% coef_dols - c(0,0,0)

inv_Sxx <- solve((t(XX)%*%XX)/nobs)
S_hat <- matrix(0,ncol=16,nrow=16)
for(j in 1:nobs)
  S_hat <- S_hat + (v_hat[j]^2 * XX[j,] %*% t(XX[j,]))/nobs

#Avar_b <- nobs * lam_hat^2 * solve((t(XX)%*%XX))
#Avar_b <- inv_Sxx %*% S_hat %*% inv_Sxx * (lam_hat/std)^2

#W_stat <- nobs * t(a_b) %*% solve( RR %*% Avar_b %*% t(RR) ) %*% a_b
#W_stat
#
#print( 1 - pchisq(W_stat,3), digits=3)

