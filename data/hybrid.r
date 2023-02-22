##  Hybrid new ISLM model
##    written by Hasumi Ryo
##      on May 17, 2007


## ------------------ 1. Defining Parameters ------------------ 

sigma <- 1.5;   # CRRA
alpha <- 0.1;
rho<-0.7;
beta <- 0.4;
q1 <- 0.4;
q2 <- 1.2;
phi <- 0.9;     # AR(1) in tech

## ------------------ Subroutines   ------------------ 
# cf. standard_rbc.r

## ------------------ 2. Steady State proc. ------------------ 
#not needed

## ------------------ 3. Model proc. ------------------ 

IS   <- expression(rho*y2a+(1-rho)*yt-sigma*((q1*ya+q2*pia+phi*vt)-pi2a)-ya); 
AS   <- expression(beta*pi2a+(1-beta)*pit+alpha*ya-pia);
dum1 <- expression(ya-y2t);
dum2 <- expression(pia-pi2t);
MPs  <- expression(phi*vt-va);        

optcon  <- list(IS, AS, dum1, dum2, MPs); 

## ------------------ 4. Solution proc. ------------------ 
xx <- c("y2a", "pi2a", "ya", "pia", "va", "y2t", "pi2t", "yt", "pit", "vt");

jopt <- rbc_jacobian(optcon, xx)

# Define Linear Coefficients  
coef <- matrix(0, nrow=length(optcon), ncol=length(xx))
for( j in 1:length(optcon))
  for (k in 1: length(xx))
    coef[j,k] <- eval(jopt[[j]][[k]]);

B <-  -coef[,1:length(jopt)];
C <-  coef[,(length(jopt)+1):(length(jopt)+length(jopt))];

A_0 <- solve(C) %*% B

A_sol <- rbc_solution(A_0)
P <- A_sol[[1]]
AA <- A_sol[[2]]

## ------------------ 5. Simulation ------------------ 

# Time span
t_sim <- 24;
# Initial Values
S1 <- matrix(c(0, 0, -0.3), ncol=1);

#SIMULATION 
X <- rbc_sim(AA, P, t_sim, S1)[[1]]

# Re-definition   
ys1 <- X[,1];
pis1 <- X[,2];
ys <- X[,1]; #Yt
pis <- X[,2]; #dp

#DRAWING FIGURES
plot(cbind(c(-3:(t_sim)) ,c(0,0,0, ys)), typ="l", xlab="Time", ylab="ys")
plot(cbind(c(-3:(t_sim)) ,c(0,0,0, pis)), typ="l", xlab="Time", col=2, ylab="pis")

