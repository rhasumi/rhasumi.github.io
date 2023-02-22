##  Basic RBC model
##    written by Hasumi Ryo
##      on May 17, 2007

## ------------------ 1. Defining Parameters ------------------ 

sigma <- 1.5   # CRRA
alpha <- 0.3   # Cob-Dag
myu <- 1       
beta <- 0.99   
delta <- 0.025 
lambda <- 2    # labor supply elasticity >1
phi <- 0.8     # AR(1) in tech


## ------------------ Subroutines   ------------------ 

# ヤコビ行列を求める関数
rbc_jacobian <- function(eq_list, x){
  # eq_list: 方程式リスト
  # x: 変数
  ans <- vector("list", length(eq_list))
  temp_list <- vector("list", length(x))
     
  for( j in 1:length(eq_list)){
    for (k in 1: length(x))
      temp_list[[k]] <- D(eq_list[[j]], x[k])
    ans[[j]] <- temp_list
  }
  return(ans)
}

# 線形化するための関数
rbc_linearise <- function(jacob, x, value){
  
  len_x <- length(x)
  len_v <- length(value)
  
  if ( identical(len_x, len_v+len_v) == F ) stop()
  
  for (j in 1: len_v){
    assign(x[j], value[j])
    assign(x[j+len_v], value[j])
  }
  
  coef <- matrix(0, nrow=length(jacob), ncol=len_x)
  for( j in 1:length(jacob))
    for (k in 1: len_x)
      coef[j,k] <- eval(jacob[[j]][[k]])
  
  TW <- matrix(value, ncol=len_v, nrow=len_v,byrow=T)
  
  B <- -coef[,1:4] * TW
  C <- coef[,5:8] * TW
  
  ans <- solve(C) %*% B 
  
  return(ans)
}

# Policy/Transition Function を求める関数
rbc_solution <- function(A){
  
  W <- eigen(A)[[2]]
  theta <- eigen(A)[[1]]
  Q <- solve(W)
  W %*% diag(theta) %*% Q
  
  # Extract stable vectors
  jw <- 0
  for (j in 1:length(theta))
    if (abs(theta[j]) >1.000000001)
      jw <- jw+1
  
  SQ <- matrix(0, ncol=length(theta), nrow=jw)
  t_jw <- 1
    for (j in 1:length(theta)){
      if (abs(theta[j]) > 1.000000001){
        SQ[t_jw,] <- Q[j,]
      	  t_jw <- t_jw + 1
      }
    }

  # Extract unstable vectors
  jjw <- 0
    for (jj in 1:length(theta))
      if (abs(theta[jj])<0.9999999999)
        jjw <- jjw+1
  
  UQ <- matrix(0, ncol=length(theta), nrow=jjw)
  t_jjw <- 1
  for (jj in 1:length(theta)){
    if (abs(theta[jj])<0.9999999999){
      UQ[t_jjw,] <- Q[jj,]
      t_jjw <- t_jjw+1
     }
  }
  
  # Extract stable roots
  VLL <- c()
  t_jjjw <- 1
  
  for (jjj in 1:length(theta)){
    if (abs(theta[jjj]) >1.0000000001){
      VLL[t_jjjw] <- theta[jjj] 
      t_jjjw <- t_jjjw+1
    }
  }
  
  # ELIMINATING UNSTABLE VECTORS
  k <- jw
  n <- jjw
  nk <- c(n, k)
  
  # Stable V (eig mat)
  VL <- solve(diag(VLL))
  
  # Elements in Q
  PA <- UQ[1:n,1:n];    PB <- UQ[1:n,(n+1):(n+k)]
  PC <- SQ[1:k,1:n];    PD <- SQ[1:k,(n+1):(n+k)]
  P <- -solve(PA) %*% PB # X(t) <- P*S(t)
  PE <- PC %*% P + PD 
  
  # SOLUTION
  PX <- solve(PE) %*% VL %*% PE
  
  return( list(P, PX, nk) )
}

# シミュレーション
rbc_sim <- function(AA, P, len_t, S_0){
  Ss <- S_0
  S <- matrix(0, nrow=len_t, ncol=ncol(AA))
  for (i in 1:len_t){
    q <- AA %*% Ss
    S[i,] <- t(q)
    Ss <- S[i,]
  }
  SY <- rbind(t(S1), S)
  X <- t(P %*% t(SY))
  return(list(X, SY))
}


## ------------------ 2. Steady State proc. ------------------ 
# SS capital & ss labor
# (1) real rate  (By SS euler)
kls <- (((1/beta)+delta-1)/alpha)^(1/(alpha-1))
# (2) wage
wstar <- (1-alpha)*(kls)^alpha
# (3) Labor and goods market clear
clstar <- kls^alpha - delta*kls
lstar <- ((wstar/myu)*(clstar^(-sigma)))^(1/(lambda+sigma))
kstar <- kls*lstar     
cstar <- clstar*lstar
vstar <- 1
Ystar <- (kstar^alpha)*(lstar^(1-alpha))


## ------------------ 3. Model proc. ------------------ 
# Eliminate Price  
# ra = (va*alpha*(ka/la)^(alpha-1)) ;
# wt = ((1-alpha)*vt*(kt/lt)^alpha);

# Optimal Conditions  & state transition
labor   <- expression(lt^lambda-((1-alpha)*vt*(kt/lt)^alpha)/(myu*ct^sigma)) # LS = LD 
euler   <- expression(ct^(-sigma) -(ca^(-sigma))*beta*(1+(va*alpha*(ka/la)^(alpha-1))-delta)) # C-Euler
capital <- expression(ka - (1-delta)*kt-vt*(kt^alpha)*(lt^(1-alpha))+ct)  # K-trans
tech    <- expression(va - phi*vt)

optcon  <- list(labor, euler, capital, tech)

# GDP (Optional)
Yt  <- expression(vt*(kt^alpha)*(lt^(1-alpha)))


## ------------------ 4. Linearization and Solution proc. ------------------ 
# Differentiation
# define vars
xx <- c("la", "ca", "ka", "va", "lt", "ct", "kt", "vt")

jopt <- rbc_jacobian(optcon, xx)

xx_star <- c(lstar, cstar, kstar, vstar)
A_0 <- rbc_linearise(jopt, xx, xx_star)

A_sol <- rbc_solution(A_0)
P <- A_sol[[1]]
AA <- A_sol[[2]]

# For GDP (optional)
xr <- c("lt", "kt", "vt")
xr_star <- c(lstar, kstar, vstar)

jy <- rbc_jacobian(Yt, xr)

for (j in 1: length(xr)){
  assign(xr[j], xr_star[j])
  assign(xr[j+length(xr)], xr_star[j])
  }

coy <- c()
for (k in 1: length(xr))
  coy[k] <- eval(jy[[1]][[k]])

ve  <- c(lstar, kstar, vstar)
NOM <- rep(Ystar, 3)
PPX <- coy * ve * NOM


## ------------------ 5. Simulation ------------------ 

# TIME&INITIAL VALUES 
# Time span

t_sim <- 24

# Initial Values
# S1 <- matrix(c(0, 0.06), ncol=1)
S1 <- matrix(c(0, 0.1), ncol=1)

# SIMULATION
temp_sim <- rbc_sim(AA, P, t_sim, S1) 
X <- temp_sim[[1]]
SY  <- temp_sim[[2]]

# Re-definition   
ci <- X[,1]
li <- X[,2]
ki <- SY[,1]
vi <- SY[,2]

XI <- matrix(c(li, ki, vi), ncol=3)
Yi <- t(PPX %*% t(XI))

# DRAWING FIGURES
plot(cbind(c(-3:(t_sim)) ,c(0,0,0, ki)), typ="l", xlab="Time", ylab="ki")
plot(cbind(c(-3:(t_sim)) ,c(0,0,0, ci)), typ="l", xlab="Time", col=2, ylab="")
lines(cbind(c(-3:(t_sim)) ,c(0,0,0, li)) ,col=4  )

