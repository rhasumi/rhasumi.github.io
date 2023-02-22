# Hayashi, Econometrics Ch. 4
# written by Hasumi Ryo
#   on July 6, 2007

greene <- read.csv("greene.csv",header=T)
g_names <- colnames(greene)      #列名の取り出し
for( j in 1:length(greene))    
  assign(g_names[j],greene[[j]]) #数値をベクトルとして列名に割り当てる

# (a)

cat(mean(Q/1000), mean(LABOR), mean(CAPITAL),mean(1-LABOR-CAPITAL),"\n")
cat(sd(Q/1000), sd(LABOR), sd(CAPITAL),sd(1-LABOR-CAPITAL),"\n")

# (b)

nobs <- length(Q)
ZZ <- cbind(rep(1,nobs), log( PL/PK ), log( PF/PK), log(Q) )

ls <- formula( LABOR ~ ZZ - 1 )
ks <- formula( CAPITAL ~ ZZ - 1 )
fs <- formula( (1-LABOR-CAPITAL) ~ ZZ - 1 )

resid_b1 <- resid( lm(ls) )
resid_b2 <- resid( lm(ks) )
resid_b3 <- resid( lm(fs) )

l_resid <- list(resid_b1, resid_b2, resid_b3)
Sig_hat <- matrix(nrow=3, ncol=3)

for(j in 1:3)
  for ( k in 1:3)
    Sig_hat[j,k] <- sum(l_resid[[j]]*l_resid[[k]])/nobs
Sig_hat
try(solve(Sig_hat))
Sig_hat[,1] + Sig_hat[,2] + Sig_hat[,3]
Sig_hat_star <- Sig_hat[c(1,3),c(1,3)]

Sig_inbv <- solve(Sig_hat_star)

sum1 <- 0
sum2 <- 0
sum3 <- 0
sum4 <- 0
sum5 <- 0
Zi <- matrix( nrow = 2, ncol = 7 )
for (j in 1:nobs){
  Zi[1,] <- c(1, 0, log(PL/PK)[j], log(PF/PK)[j],            0 , log(Q)[j],        0 )
  Zi[2,] <- c(0, 1,            0 , log(PL/PK)[j], log(PF/PK)[j],         0, log(Q)[j])
  yi <- rbind(LABOR[j],(1-LABOR-CAPITAL)[j])
  xi <- c(1, log(PL/PK)[j], log(PF/PK)[j], log(Q)[j])
  sum1 <- sum1 + (t(Zi) %*% Sig_inv %*% Zi)/nobs
  sum2 <- sum2 + (t(Zi) %*% Sig_inv %*% yi)/nobs
  sum3 <- sum3 + (xi %*% t(xi))/nobs
  sum4 <- sum4 + (yi %x% xi)/nobs
  sum5 <- sum5 + (Zi %x% xi)/nobs
}

avar_d <- solve(sum1)
d_re <- avar_d %*% ( sum2 )
d_re
sqrt(diag(avar_d/nobs))

table44 <- matrix(ncol=3, nrow=12)

rownames(table44) <- c("a1", "a2", "a3", "g11", "g12", "g13", "g22", "g23", "g33", "g1q", "g2q", "g3q")
table44[c(1,3,4,6,9,10,12),1] <- d_re
table44[c(1,3,4,6,9,10,12),2] <- sqrt(diag(avar_d/nobs))

table44[2,1] <- 1 - table44[1,1] - table44[3,1] 
  table44[2,2] <- sqrt((t(c(1,1))%*%avar_d[c(1,2),c(1,2)]%*%c(1,1))/nobs)
table44[5,1] <- - table44[4,1] - table44[6,1]
  table44[5,2] <- sqrt((t(c(1,1))%*%avar_d[c(3,4),c(3,4)]%*%c(1,1))/nobs)
table44[8,1] <- - table44[6,1] - table44[9,1]
  table44[8,2] <- sqrt((t(c(1,1))%*%avar_d[c(4,5),c(4,5)]%*%c(1,1))/nobs)
table44[7,1] <- - table44[5,1] - table44[8,1]
  table44[7,2] <- sqrt((t(c(1,2,1))%*%avar_d[c(3,4,5),c(3,4,5)]%*%c(1,2,1))/nobs)
table44[11,1] <- - table44[10,1] - table44[12,1];
  table44[11,2] <- sqrt((t(c(1,1))%*%avar_d[c(6,7),c(6,7)]%*%c(1,1))/nobs)

table44[,3] <- table44[,1]/table44[,2]
print(table44, digits=3)

# (c)

g_d_re <- (sum4) - (sum5) %*% d_re
S_d_re <- nobs * ( t(g_d_re) %*% solve( Sig_hat_star %x% sum3 ) %*% g_d_re )
1-pchisq(S_d_re, 2*4-7)

# (e)

RR <- c(0,0,1,0,0,-1,0,0)
dd <- c(coef(lm(ls)),coef(lm(fs)))

avar_nr <- Sig_hat_star %x% solve(t(ZZ) %*% ZZ / nobs)
W_stat <- nobs * (dd[3] - dd[6]) %*% solve( RR %*% avar_nr %*% RR ) %*%  (dd[3] - dd[6])

print(c(S_d_re, W_stat))

# (f)

fL <- ZZ %*% d_re[c(1,3,4,6)]
fF <- ZZ %*% d_re[c(2,4,5,7)]
fK <- 1 - fL - fF
mean(table44[5,1]/(fL*fK) + 1)

