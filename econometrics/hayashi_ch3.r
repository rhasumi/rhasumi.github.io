# Hayashi, Econometrics Ch. 3
# written by Hasumi Ryo
#   on July 4, 2007

grilic <- read.csv("grilic.csv",header=T)
names <- colnames(grilic)      #列名の取り出し
for( j in 1:length(grilic))    
  assign(names[j],grilic[[j]]) #数値をベクトルとして列名に割り当てる

# (a)

mean(grilic)
sd(grilic)
cor(IQ,S)

# (b)
library(systemfit) # 予めダウンロードする必要があります

nobs <- length(grilic[[1]])

dumyear <- matrix(0, ncol=7, nrow=nobs)
dumyear[which(YEAR==66),1] <- 1
dumyear[which(YEAR==67),2] <- 1
dumyear[which(YEAR==68),3] <- 1
dumyear[which(YEAR==69),4] <- 1
dumyear[which(YEAR==70),5] <- 1
dumyear[which(YEAR==71),6] <- 1
dumyear[which(YEAR==73),7] <- 1

h <- cbind(EXPR, TENURE, RNS, SMSA, dumyear)

line_1 <- lm(LW ~ S + h - 1 )
line_2 <- lm(LW ~ S + IQ + h - 1 )

summary(line_1)
summary(line_2)

excl <- cbind(MED,KWW,MRT,AGE)
inst_3 <- IQ ~ S + h + excl - 1
eqLW <- LW ~ S + IQ + h - 1

line_3 <- systemfit( "2SLS", list(eqLW), inst = inst_3 )
summary(line_3)

# (c)

res <- resid( line_3 )[[1]]
XX <- cbind( S, h, excl )

Sargan <- function(X,eps,n){
  ans <- ( t(eps) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% eps ) / (t(eps) %*% eps) * n
  return( ans )
}
S_stat <- Sargan( XX, res, nobs )
S_stat
1-pchisq(S_stat, 3)

# (d)

IQ_hat <- predict(lm(inst_3))
eq_d <-  LW ~ S + IQ_hat + h - 1

result_d <- lm(eq_d)
summary(result_d)

# (e)

inst_4a <- IQ ~ + h + excl - 1
inst_4b <- S ~ + h + excl - 1
eqLW_4 <- LW ~ S + IQ + h - 1

line_4 <- systemfit( "2SLS", list(eqLW), inst = list(inst_4a,inst_4b) )
summary(line_4)

e_hat <- resid( line_4 )[[1]]
xx <- cbind( h, excl )

S_stat_f <- Sargan( xx, e_hat, nobs )
S_stat_f
print(1-pchisq(S_stat_f, 2),digits=20)


# GMM (line 5)

zz <- cbind(S, IQ, h)
yy <- LW

S_hat <- matrix(0, ncol=length(xx[1,]), nrow=length(xx[1,]))
for (j in 1:nobs)
  S_hat <- S_hat + e_hat[j]^2 * xx[j,] %*% t(xx[j,])/nobs

dim_xz <- dim(xx[1,] %*% t(zz[1,]))

S_xz <- matrix(0, nrow=dim_xz[1], ncol=dim_xz[2] )
for (j in 1:nobs)
  S_xz <- S_xz +  xx[j,] %*% t(zz[j,])/nobs


s_xy <- matrix(0, nrow=length(xx[1,]), ncol=1 )
for (j in 1:nobs)
  s_xy <- s_xy +  xx[j,] * yy[j]/nobs

W_hat <- solve(S_hat)

d_gmm <- solve( t(S_xz) %*% W_hat %*% S_xz ) %*% t(S_xz) %*% W_hat %*% s_xy
d_gmm
Avar_d <- solve(t(S_xz)%*%W_hat%*%S_xz)
sqrt( diag(Avar_d)/nobs )
sqrt(sum((yy - zz %*% d_gmm)^2)/nobs)

g_n <- s_xy - S_xz %*% d_gmm

J_stat <- nobs * t(g_n) %*% W_hat %*% g_n
J_stat
print(1-pchisq(J_stat, 2),digits=20)


# (f)

xx2 <- cbind(h, excl,S)
zz <- cbind(S, IQ, h)
yy <- LW

e_hat2 <- resid( line_3 )[[1]]
S_hat2 <- matrix(0, ncol=length(xx2[1,]), nrow=length(xx2[1,]))
for (j in 1:nobs)
  S_hat2 <- S_hat2 + e_hat2[j]^2 * xx2[j,] %*% t(xx2[j,])/nobs

dim_xz2 <- dim(xx2[1,] %*% t(zz[1,]))

S_xz2 <- matrix(0, nrow=dim_xz2[1], ncol=dim_xz2[2] )
for (j in 1:nobs)
  S_xz2 <- S_xz2 +  xx2[j,] %*% t(zz[j,])/nobs

s_xy2 <- matrix(0, nrow=length(xx2[1,]), ncol=1 )
for (j in 1:nobs)
  s_xy2 <- s_xy2 +  xx2[j,] * yy[j]/nobs

W_hat2 <- solve(S_hat2)

d_gmm2 <- solve( t(S_xz2) %*% W_hat2 %*% S_xz2 ) %*% t(S_xz2) %*% W_hat2 %*% s_xy2
d_gmm2

g_n2 <- s_xy2 - S_xz2 %*% d_gmm2

J_stat2 <- nobs * t(g_n2) %*% W_hat2 %*% g_n2
J_stat2

W_hat1 <- solve(S_hat2[1:15,1:15])
d_gmm1 <- solve( t(S_xz2[1:15,]) %*% W_hat1 %*% S_xz2[1:15,] ) %*% t(S_xz2[1:15,]) %*% W_hat1 %*% s_xy2[1:15,]

g_n1 <- s_xy2[1:15] - S_xz2[1:15,] %*% d_gmm1


J_stat1 <- nobs * t(g_n1) %*% W_hat1 %*% g_n1
J_stat1

C_stat <- J_stat2 - J_stat1
C_stat


# (g)

inst_ga <- IQ ~ + h + MRT + AGE - 1
inst_gb <- S ~ + h + MRT + AGE - 1
eqLW_g <- LW ~ S + IQ + h - 1

line_g <- systemfit( "2SLS", list(eqLW), inst = list(inst_ga,inst_gb) )
summary(line_g)

