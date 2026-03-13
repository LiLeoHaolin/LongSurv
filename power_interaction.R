power <- function(theta, beta1, beta2, alpha.mean, alpha.var, n, event.prop.null, EA, p, ui, sigma.square, trans, sig){
  # moment generating functions
  mom1 <- function(mean, var, theta){
    x <- exp(t(mean) %*% theta + 0.5*t(theta) %*% var %*% theta)
    return(x)
  }
  mom2 <- function(mean, var, theta){
    x <- (mean + var %*% theta) * as.numeric(exp(t(mean) %*% theta + 0.5*t(theta) %*% var %*% theta))
    return(x)
  }
  mom3 <- function(mean, var, theta){
    x <- as.numeric(exp(t(mean) %*% theta + 0.5*t(theta) %*% var %*% theta))*var + as.numeric(exp(t(mean) %*% theta + 0.5*t(theta) %*% var %*% theta))*(mean + var %*% theta) %*% t(mean + var%*%theta)
    return(x)
  }
  
  moment1_null <- mom1(alpha.mean, alpha.var, theta)
  moment2_null <- mom2(alpha.mean, alpha.var, theta)
  moment3_null <- mom3(alpha.mean, alpha.var, theta)
  
  moment1_alt <- as.numeric(exp(beta1))*mom1(alpha.mean, alpha.var, theta+beta2)
  moment2_alt <- as.numeric(exp(beta1))*mom2(alpha.mean, alpha.var, theta+beta2)
  moment3_alt <- as.numeric(exp(beta1))*mom3(alpha.mean, alpha.var, theta+beta2)
  
  # prereq quantities
  Ui <- numeric()
  for (i in 1:p){
    Ui <- cbind(Ui, ui^(i-1))
  }
  event.prop <- as.numeric((1-EA)*event.prop.null + EA*moment1_alt/moment1_null*event.prop.null)
  
  # sigma.one component
  s0 <- EA*moment1_alt+(1-EA)*moment1_null
  s1.div.s0 <- rbind((EA*moment2_alt+(1-EA)*moment2_null), EA*moment1_alt, EA*moment2_alt)/as.numeric(s0)
  s2 <- cbind(rbind((EA*moment3_alt+(1-EA)*moment3_null), EA*t(moment2_alt), EA*moment3_alt), rbind(EA*moment2_alt, EA*moment1_alt, EA*moment2_alt), rbind(EA*moment3_alt, EA*t(moment2_alt), EA*moment3_alt))
  sigma.one <- (s2/as.numeric(s0) - s1.div.s0 %*% t(s1.div.s0))*event.prop
  
  # H1 component
  H1 <- cbind(rbind(sigma.square*solve(t(Ui) %*% Ui), matrix(0, nrow=1, ncol=p), EA*sigma.square*solve(t(Ui) %*% Ui)), matrix(0, nrow=2*p+1, ncol=1), rbind(EA*sigma.square*solve(t(Ui) %*% Ui), matrix(0, nrow=1, ncol=p), EA*sigma.square*solve(t(Ui) %*% Ui)))*event.prop
  
  # Asymptotic variance
  H = H1 
  #H = H1 + H2 + H3 - H6
  Sigma <- solve(sigma.one) + solve(sigma.one) %*% H %*% solve(sigma.one) 
  
  # NCP
  xi <- rbind(theta, beta1, beta2)
  NCP <- as.numeric(t(trans %*% xi) %*% solve(trans %*% (Sigma/n) %*% t(trans)) %*% (trans %*% xi))
  #NCP <- as.numeric(t(trans %*% xi) %*% solve(trans %*% (Sigma) %*% t(trans)) %*% (trans %*% xi))
  
  # power
  power <- 1-pchisq(qchisq(1-sig, p), p, NCP)
  
  # output
  return(power)
  #return(event.prop)
}

power(theta = matrix(0.25), # longitudinal effect
      beta1 = matrix(-0.05), # trt effect under alternative
      beta2 = matrix(-0.15), # interaction under alternative
      alpha.mean = matrix(6.7), # mean of random effect
      alpha.var = matrix(20.25), # variance of random effect
      n = 4229, # sample size
      event.prop = 0.054, # event proportion
      #event.prop = 0.135014,
      EA = 0.5, # treatment allocation
      p = 1, # number of random effect
      ui <- matrix(c(-1), nrow = 1), # time measures for longitudinal biomarkers
      sigma.square = 1, # variance of random error
      trans = matrix(c(0,0,1), nrow = 1), # transfomation matrix of parameters for hypothesis testing
      sig = 0.05) # significance level