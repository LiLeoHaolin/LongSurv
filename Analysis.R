library(Rcpp)
library(RcppArmadillo)
setwd("/nas/longleaf/home/haolin/long_surv/rcpp")
sourceCpp('RcppFunctions.cpp')

######################################
# start the simulation
######################################

#path <- "C:/Users/haolin/Desktop/p01/"
#path <- "C:/Users/haoli/Desktop/LongSurv/05_simulated_data/"
path <- "/nas/longleaf/home/haolin/long_surv/data_v2/05_simulated_data/"
#p <- 3
p <- 1
sigma.hat <- 1
rej.cum = 0
nsim = 1000

beta1 <- 0
beta2 <- matrix(0)
theta <- matrix(-0.2)

H = matrix(c(0,0,1), nrow = 1)

for (z in 1:nsim){
  # z=2
  cat(z)
  
  try({

dat <- read.csv(file = paste0(path, 'dat', z,'.csv'))
U <- read.csv(file = paste0(path, 'U', z,'.csv'))
U <- as.matrix(U)
dat$m <- (dat$X1!=0)+(dat$X2!=0)+(dat$X3!=0)+(dat$X4!=0)+(dat$X5!=0)+(dat$X6!=0)+(dat$X7!=0)+(dat$X8!=0)+(dat$X9!=0)+(dat$X10!=0)
m <- dat$m
max.m <- max(m)
    
######################################
# estimation
######################################

# solve for alpha
alpha.hat <- matrix(0, nrow = nrow(dat), ncol = p, byrow = T)
for (i in 1:nrow(dat)){
  xi <- t(dat[i, -c(1, 2, 3, 4)])
  Ui <- as.matrix(U[c(((i-1)*max.m+1):(i*max.m)), -1])
  xi <- xi[1:m[i],]
  Ui <- Ui[1:m[i],]
  if (p==1){
    linearMod <- lm(xi ~ 1)
  } else{
    linearMod <- lm(xi ~ Ui)
  }
  #linearMod <- lm(xi ~ 1) # random intercept only model
  #linearMod <- lm(xi ~ Ui) # intercept already included
  alpha.hat[i, ] <- linearMod$coefficients
}
    
    # solve for sigma
    # dat.sigma <- data.frame(id=rep(1:N, each = max.m), out=matrix(t(X), nrow = N*max.m))
    # dat.sigma <- cbind(dat.sigma, U)
    # dat.sigma <- dat.sigma[dat.sigma$out!=0, ]
    # print(nrow(dat.sigma))
    # sum(m)
    ### need some human modifications to the specification to the model.
    # rand.eff.mod <- lmer(out ~  (`2`+`3`|id), data = dat.sigma)
    # sigma.hat <- summary(rand.eff.mod)$sigma
    
trt <- as.matrix(dat$trt)
time <- as.matrix(dat$time)
status <- as.matrix(dat$status)
n <- nrow(dat)

# corrected score
corr.score = function(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p){
  U.new = U
  trt.new = trt
  time.new = time
  status.new = status
  theta.new = theta
  beta1.new = beta1
  beta2.new = beta2
  sigma.new = sigma
  t.new = t
  alpha.new = alpha
  comp1 <- 0
  comp2 <- 0
  comp3 <- 0
  for (j in 1:n) comp1 <- comp1 +1*(status.new[j] == 1)*(time.new[j]<=t.new) * (t(t(alpha.new[j,])) 
       - S_tilde_oneone_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[j], alpha = alpha.hat, maxm = max.m, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[j], alpha = alpha.hat, maxm = max.m, p = p)))
  for (k in 1:n) comp2 <- comp2 +1*(status.new[k] == 1)*(time.new[k]<=t.new) * (trt.new[k] 
       - S_tilde_onetwo_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[k], alpha = alpha.hat, maxm = max.m, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[k], alpha = alpha.hat, maxm = max.m, p = p)))
  for (l in 1:n) comp3 <- comp3 +1*(status.new[l] == 1)*(time.new[l]<=t.new) * (t(t(alpha.new[l,]))*as.numeric(trt.new[l]) 
       - S_tilde_onethree_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[l], alpha = alpha.hat, maxm = max.m, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=time.new[l], alpha = alpha.hat, maxm = max.m, p = p)))
  return(rbind(comp1, as.matrix(comp2), comp3))
}
# test
# corr.score(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma.hat, t=100, alpha = alpha.hat, maxm = max.m, p = p)

# corr.score(U=U, dat=dat, beta1=beta1, theta=theta, beta2=beta2, sigma = sigma.hat, t=100, alpha = alpha.hat)


#######################
### estimation way1: NR
    
tol = 10^-4
Theta.hat <- as.matrix(rbind(theta, beta1, beta2))
#Theta.hat[(p+1)] <- 0.5
#Theta.hat <- matrix(0, nrow = (2*p+1))
maxit = 50
iter = 0
eps = Inf
start = Sys.time()
while(eps > tol & iter < maxit){
  # save the previous value
  Theta.hat.0 = Theta.hat
  theta.hat <- as.matrix(Theta.hat[0:p])
  beta1.hat <- as.numeric(Theta.hat[p+1])
  beta2.hat <- as.matrix(Theta.hat[(p+2):(2*p+1)])
  # calculate h, the increment
  Score <- corr.score(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=max(time), alpha = alpha.hat, maxm = max.m, p = p)
  Sigma1.hat <- 0
  for (i in 1:nrow(dat)){
    Sigma1.hat <- Sigma1.hat + (1/n)* (status[i] == 1)* 
      ((rbind(S_tilde_one_prime_one_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
              S_tilde_one_prime_two_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
              S_tilde_one_prime_three_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)))
       /as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p))
       - rbind(S_tilde_oneone_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
               S_tilde_onetwo_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
               S_tilde_onethree_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)) 
       %*% t(S_tilde_zero_prime_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)) 
       /((as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)))^2))
  }
  Info <- -n*Sigma1.hat
  h = - solve(Info) %*% Score
  # update beta
  Theta.hat = Theta.hat + h
  # calculate the diff in beta, could also use the log likelihood if we wanted
  eps  = t(Theta.hat - Theta.hat.0) %*% (Theta.hat - Theta.hat.0)
  # update the iteration number
  iter = iter + 1
  if(iter == maxit) warning("Iteration limit reached without convergence")
  # print out info to keep track
  cat(sprintf("Iter: %d Theta: %f h: %f eps:%f\n",iter,Theta.hat, h, eps))
}
    
# format theta, beta1, and beta2
theta.hat <- as.matrix(Theta.hat[0:p])
beta1.hat <- as.numeric(Theta.hat[p+1])
beta2.hat <- as.matrix(Theta.hat[(p+2):(2*p+1)])
  
    
    

######################################
# var est
######################################

# definition of E.tilde
E.tilde <- function(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p){
  comp1 <- S_tilde_oneone_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p))
  comp2 <- S_tilde_onetwo_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p))
  comp3 <- S_tilde_onethree_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p)/as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta, beta1 = beta1, beta2=beta2, sigma = sigma, t=t, alpha = alpha, maxm = maxm, p = p))
  return(rbind(comp1, as.matrix(comp2), comp3))
}
# test
# E.tilde(U=U, trt=trt, time=time, status=status, n=n, beta1=beta1.hat, theta=theta.hat, beta2=beta2.hat, sigma = sigma.hat, t=0.5, alpha = alpha.hat, maxm = max.m, p = p)

# definition of xi
xi <- function(i){
  Ui <- as.matrix(U[c(((i-1)*max.m+1):(i*max.m)), ])
  Ui <- Ui[1:m[i],]
  Di <- t(theta.hat+beta2.hat*dat$trt[i]) %*% solve(t(Ui)%*%Ui) %*% (theta.hat+beta2.hat*dat$trt[i])
  sum.2 <- 0
  for (j in 1:nrow(dat)){
    sum.2 <- sum.2 + as.numeric((1/nrow(dat))*dat$status[j]* (dat$time[i] >= dat$time[j])* exp(beta1.hat*dat$trt[i]) * exp(t(theta.hat+beta2.hat*dat$trt[i]) %*% t(t(alpha.hat[i,]))-0.5*sigma.hat^2*Di))/
    as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[j], alpha = alpha.hat, maxm = max.m, p = p)) * (rbind((t(t(alpha.hat[i,]))-(sigma.hat^2)*solve(t(Ui)%*%Ui)%*%(theta.hat + beta2.hat*dat$trt[i])), as.matrix(dat$trt[i]), dat$trt[i] * ((t(t(alpha.hat[i,]))-(sigma.hat^2)*solve(t(Ui)%*%Ui)%*%(theta.hat + beta2.hat*dat$trt[i])))) 
       - E.tilde(U=U, trt=trt, time=time, status=status, n=n, beta1=beta1.hat, theta=theta.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[j], alpha = alpha.hat, maxm = max.m, p = p))
  }
  first <- 1*(dat$status[i] == 1)* (rbind(t(t(alpha.hat[i,])), as.matrix(dat$trt[i]), dat$trt[i]*t(t(alpha.hat[i,]))) 
                                    - E.tilde(U=U, trt=trt, time=time, status=status, n=n, beta1=beta1.hat, theta=theta.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p))-sum.2
  return(first)
}
# test
# xi(3)%*%t(xi(3))

# Sigma2.hat
Sigma2.hat <- matrix(0, nrow=2*p+1, ncol=2*p+1)
for (i in 1:nrow(dat)){
  cat(i)
  Sigma2.hat <- Sigma2.hat+ xi(i)%*%t(xi(i))
}
Sigma2.hat <- (1/nrow(dat))*Sigma2.hat

# Sigma1.hat
Sigma1.hat <- 0
for (i in 1:nrow(dat)){
  cat(i)
  Sigma1.hat <- Sigma1.hat + (1/nrow(dat))* (dat$status[i] == 1)* ((rbind(
    S_tilde_one_prime_one_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p),
    S_tilde_one_prime_two_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
    S_tilde_one_prime_three_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)))/
      as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p))- rbind(
        S_tilde_oneone_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
        S_tilde_onetwo_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p), 
        S_tilde_onethree_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)) %*% t(S_tilde_zero_prime_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)) /(as.numeric(S_tilde_zero_Rcpp(U = U, trt = trt, time = time, status = status, n = n, theta=theta.hat, beta1 = beta1.hat, beta2=beta2.hat, sigma = sigma.hat, t=time[i], alpha = alpha.hat, maxm = max.m, p = p)))^2)
}

write.csv(Sigma1.hat,file=paste0('/nas/longleaf/home/haolin/long_surv/round2/results/matrix1_', z, '.csv'))
Sigma.hat <- solve(Sigma1.hat) %*% Sigma2.hat %*% solve(Sigma1.hat)/nrow(dat)
write.csv(Sigma2.hat-Sigma1.hat,file=paste0('/nas/longleaf/home/haolin/long_surv/round2/results/matrix', z, '.csv'))

######################################
# Summary
######################################

Var.hat <- diag(Sigma.hat)
se <- sqrt(Var.hat)

result <- data.frame(var.n=1:(2*p+1), true = rbind(theta, beta1, beta2), est = Theta.hat, se = se)
result$inc <- 1*(result$true <= result$est+1.96*result$se)*(result$true >= result$est-1.96*result$se)
result$nsim <- z
write.csv(result,file=paste0('/nas/longleaf/home/haolin/long_surv/round2/results/result', z, '.csv'))

Sigma.hat <- as.matrix(Sigma.hat)
est <- as.matrix(result$est)

W = as.numeric(t(H %*% est) %*% solve(H %*% Sigma.hat %*% t(H)) %*% (H %*% est))
rej = as.numeric(W>qchisq(0.95, 1))

output <- data.frame(W=W, rej=rej)
#write.csv(result,file=paste0('C:/Users/haoli/Desktop/LongSurv/01_simulation/results/testing/result', z, '.csv'))
write.csv(output,file=paste0('/nas/longleaf/home/haolin/long_surv/round2/results/output', z, '.csv'))

  })
}






