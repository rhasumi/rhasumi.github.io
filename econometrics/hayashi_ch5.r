#
# Hayashi Ch. 5
#

# Empirical Exercise

sum_hes <- read.table("sum_hes.asc",header=F,fill=T)

# (a)

y <- log(sum_hes[[3]][27*c(1:125)])
X <- matrix(1,ncol=6,nrow=125)

y_0 <- log(sum_hes[[3]][2 + 27*c(0:124)])
X[,2] <- y_0

s <- c()
for(j in 1:125)
  s[j] <- mean( sum_hes[[4]][c(2:27)+27*(j-1)] )
X[,3] <- log(s)

n <- c()
for(j in 1:125)
 n[j] <- mean(sum_hes[[2]][c(3:27)+27*(j-1)]/sum_hes[[2]][c(2:26)+27*(j-1)]-1)*100
X[,4] <- log(n+5)

X[,5] <- sum_hes[[2]][1 + 27*c(0:124)]
X[,6] <- sum_hes[[3]][1 + 27*c(0:124)]

plot(cbind(y_0,y-y_0))

ols <- lm(y ~ X -1)
summary(ols)

eps <- resid(ols)
nobs <- 125

Sxx <- matrix(0,nrow=6,ncol=6)
S <- matrix(0,nrow=6,ncol=6)

for(j in 1:nobs){
  Sxx <- Sxx + X[j,] %*% t(X[j,]) / nobs 
  S <- S + eps[j]*X[j,] %*% t(eps[j]*X[j,]) / nobs
}

Avar <- solve(Sxx) %*% S %*% solve(Sxx)
Avar_rho <- sqrt(Avar[2,2]/nobs)
rho_hat <- coef(ols)[2]

lambda_hat <- - log(rho_hat)/25
Avar_lambda <- Avar_rho / (25*rho_hat)^2

cat(lambda_hat,Avar_lambda,"\n")

# (b)

XX <- matrix(0, ncol=25+1,nrow=125*25)
for(j in 1:25)
  XX[c(0:124)*25+j,j] <- 1

yy <- c()
for(jj in 1:125){
  y_l <- log(sum_hes[[3]][c(3:27)+(jj-1)*27])
  y_r <- log(sum_hes[[3]][c(2:26)+(jj-1)*27])
  yy[c(1:25)+(jj-1)*25] <- y_l - mean(y_l)
  XX[c(1:25)+(jj-1)*25,26] <- y_r - mean(y_r)
}

ols_b <- lm(yy~XX -1)
summary(ols_b)

rho_hat_b <- coef(ols_b)[26]
-log(rho_hat_b)


