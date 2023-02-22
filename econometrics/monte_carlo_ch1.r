# Monte Carlo Simulation
# written by Hasumi Ryo
#   on June 15, 2007

beta_1 <- 1
beta_2 <- 0.5
sigs <- 1

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

# simulation A (ristricted)
nrep_a <- 10^3
c_a <- 0
sim_b1_a <- c()
sim_b2_a <- c()

x_a <-  x()

for( j in 1:nrep_a){
  y <- beta_1 + x_a*beta_2 + rnorm(n) / sqrt(sigs)
  X <- cbind( rep(1,n), x_a)
  olse <- lm(y ~ x_a)
  ss <- deviance(olse)
  sim_b1_a[j] <- olse[[1]][1]
  sim_b2_a[j] <- olse[[1]][2]
  if ( abs(sim_b2_a[j]-beta_2)/ summary(olse)[[4]][2,2] < abs(qt(0.025, n-2)) )
    c_a <- c_a + 1
}

# print
c_a
c_a/nrep_a
mean(sim_b1_a)
mean(sim_b2_a)

# simulation B (not restricted)

nrep_b <- 10^3
c_b <- 0
sim_b1_b <- c()
sim_b2_b <- c()

for( j in 1:nrep_b){
  x_b <- x()
  y <- beta_1 + x_b*beta_2 + rnorm(n) / sqrt(sigs)
  olse <- lm(y ~ x_b)
  ss <- deviance(olse)
  sim_b1_b[j] <- olse[[1]][1]
  sim_b2_b[j] <- olse[[1]][2]
  if ( abs(sim_b2_b[j]-beta_2)/ summary(olse)[[4]][2,2] < abs(qt(0.025, n-2)) )
    c_b <- c_b + 1
}

# print
c_b
c_b/nrep_b
mean(sim_b1_b)
mean(sim_b2_b)

