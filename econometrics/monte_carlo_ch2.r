# Monte Carlo Simulation
# written by Hasumi Ryo
#   on May 8, 2007

# simulation 1

beta_1 <- 1
beta_2 <- 0.5

n <- 32
cc <- 2
phi <- 0.6

r <- cumprod(rep(phi,n))/phi

A <- matrix(0, ncol=n, nrow=n)
for( j in 1:n){
  A[j:n,j] <- cumprod(rep(phi,n+1-j))/phi
}

d <- c()
for( j in 1:n){
d[j] <- sum(cumprod(rep(phi,j))/phi)*cc
}

x <- function(){
  x_0 <- rnorm(1)/sqrt(1-phi^2) + cc
  ans <- r*x_0 + d + A%*% rnorm(n)
  return(ans)
}


nrep <- 10^4
c_1 <- 0
c_2 <- 0
sim_b1 <- c()
sim_b2 <- c()

for( j in 1:nrep){
  xx <- x()
  y <- beta_1 + xx*beta_2 + runif(n) - 0.5
  olse <- lm(y ~ xx)
  ss <-  deviance(olse)
  sim_b1[j] <- olse[[1]][1]
  sim_b2[j] <- olse[[1]][2]
  if ( abs(sim_b2[j]-beta_2)/ summary(olse)[[4]][2,2] < abs(qt(0.025, 30)) )
    c_1 <- c_1 + 1
  if ( abs(sim_b2[j]-beta_2)/ summary(olse)[[4]][2,2] < abs(qnorm(0.025) ) )
    c_2 <- c_2 + 1
}

# print
cat(c_1,c_2,"\n")
cat(c_1/nrep,c_2/nrep,"\n")
mean(sim_b1)
mean(sim_b2)


# simulation 2
n <- 50
p <- 4
nrep <- 10^4

c_bp <- 0
c_lb <- 0

for( j in 1:nrep){
  z <- runif(50) - 0.5
  gam <- c()
  for (j in 1:p)
    gam[j] <- sum( (z[1:(n-j)]-mean(z))*(z[(1+j):n]-mean(z)) )/n
  #var_z <- var(z)
  var_z <- sum((z[1:n]-mean(z))^2)/n
  bp_Q <- n*sum( (gam/var_z)^2 )
  lb_Q <- n*(n+2) * sum( (gam/var_z)^2 / (n-c(1:p)) )
  if ( bp_Q < qchisq(0.95,p) )
    c_bp <- c_bp + 1
  if ( lb_Q < qchisq(0.95,p) )
    c_lb <- c_lb + 1
}

# print
cat(c_bp,c_lb,"\n")
cat(c_bp/nrep,c_lb/nrep,"\n")

