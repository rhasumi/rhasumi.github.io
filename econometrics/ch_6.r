#
# Hayashi Ch. 6
#

dm <- read.table("dm.asc", header=F)
pound <- read.table("pound.asc", header=F)
yen <- read.table("yen.asc", header=F)

dm_s <- log(dm[[2]])
dm_f <- log(dm[[3]])
dm_s30 <- log(dm[[4]])

pound_s <- log(pound[[2]])
pound_f <- log(pound[[3]])
pound_s30 <- log(pound[[4]])

yen_s <- log(yen[[2]])
yen_f <- log(yen[[3]])
yen_s30 <- log(yen[[4]])

# (1b)

dm_e <- dm_s30 - dm_f
acf(dm_e, lag.max=40)

pound_e <- pound_s30 - pound_f
acf(pound_e, lag.max=40)

yen_e <- yen_s30 - yen_f
acf(yen_e, lag.max=40)

# (1c)

dm_ds <- (dm_s[-1]-dm_s[-length(dm_s)])
acf(dm_ds, lag.max=40)
Box.test(dm_ds, lag = 40, type = c("Ljung-Box"))

pound_ds <- (pound_s[-1]-pound_s[-length(pound_s)])
acf(pound_ds, lag.max=40)
Box.test(pound_ds, lag = 40, type = c("Ljung-Box"))

yen_ds <- (yen_s[-1]-yen_s[-length(yen_s)])
acf(yen_ds, lag.max=40)
Box.test(yen_ds, lag = 40, type = c("Ljung-Box"))

# (1d) for yen/$

nobs <- length(yen_s)

mean( yen_s30 - yen_s )*1200; sd( yen_s30 - yen_s)*1200
mean( yen_f - yen_s )*1200; sd( yen_f - yen_s )*1200
mean( yen_e )*1200; sd(yen_e)*1200


m_ye <- mean(yen_e)

gamma_0 <- sum((yen_e-m_ye)^2)/nobs
gamma_1 <- sum( (yen_e-m_ye)[1:(nobs-1)] * (yen_e-m_ye)[2:nobs])/nobs
gamma_2 <- sum( (yen_e-m_ye)[1:(nobs-2)] * (yen_e-m_ye)[3:nobs])/nobs
gamma_3 <- sum( (yen_e-m_ye)[1:(nobs-3)] * (yen_e-m_ye)[4:nobs])/nobs
gamma_4 <- sum( (yen_e-m_ye)[1:(nobs-4)] * (yen_e-m_ye)[5:nobs])/nobs

sqrt((gamma_0 + 2*(gamma_1+gamma_2+gamma_3+gamma_4))/nobs)*1200

# (1f)

yy <- (yen_s30 - yen_s)*1200
xx <- cbind(rep(1,nobs),(yen_f-yen_s)*1200)
ols_yen <- lm(yy~xx-1)
summary(ols_yen)

Sxx <- matrix(0, ncol=2, nrow=2)
for( j in 1:nobs)
  Sxx <- Sxx + xx[j,] %*% t(xx[j,]) / nobs

eps <- as.double(resid( ols_yen ))

g <- xx
g[,1] <-  xx[,1]*eps
g[,2] <-  xx[,2]*eps

qq <- 12

gamma <- vector("list", qq+1)
for( j in 0:qq ){
   temp <- matrix(0, ncol=2, nrow=2)
   for (tt in 1:(nobs-j) )
     temp <- temp +  g[tt,] %*% t(g[tt+j,]) /nobs
     gamma[[j+1]] <- temp
}

SS <- gamma[[1]]

for (j in 1:qq){
  # print( 1 - j/(qq+1) )
  SS <- SS + (1 - j/(qq+1)) * (gamma[[j+1]]+t(gamma[[j+1]]))
}

Avar_b <- solve(Sxx) %*% SS %*% solve(Sxx)
sqrt(Avar_b[1,1]/nobs)
sqrt(Avar_b[2,2]/nobs)

RR <-  matrix(c(1,0,0,1),ncol=2)
a_b <- RR %*% coef(ols_yen) - c(0,1)

W_stat <- nobs * t(a_b) %*% solve( RR %*% Avar_b %*% t(RR) ) %*% a_b
W_stat

print( 1 - pchisq(W_stat,2), digits=10)


# (2b)

mishkin <- read.table("MISHKIN.asc",header=F)
len <- length(mishkin[[1]])
pai3 <- mishkin[[4]]; p_t3 <- ts( pai3, start=c(1950,2), frequency=12 )
tb3  <- mishkin[[6]]; R_t <- ts( tb3, start=c(1950,2), frequency=12 )
r_t3 <-  R_t - p_t3


r_t3_c <- window(r_t3,c(1953,1),c(1971,7))
nobs <- length(r_t3_c)
nobs

acf(as.double(r_t3_c))

# (2d)

yy <- window(p_t3,c(1953,1),c(1971,7))
xx <- cbind(rep(1,nobs), window(R_t,c(1953,1),c(1971,7)))

ols_2d <- lm(yy ~ xx - 1)
summary( ols_2d )


