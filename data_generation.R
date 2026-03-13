######################################
# start the simulation
######################################

#path <- "C:/Users/haoli/Desktop/LongSurv/01_simulation/naive_cox/data/"
path <- "/nas/longleaf/home/haolin/long_surv/data_v3/05_simulated_data/"
#path <- "C:/Users/haoli/Desktop/LongSurv/05_simulated_data/"

nsim = 1000

for (z in 1:nsim){
  # z=1
  cat(z)
  try({
    
    ######################################
    # quantities related to the set up
    ######################################
    
    # N = sample size
    N <- 1500
    
    # dimension
    # p <- 3
    p <- 1
    
    # parameters related to censoring time
    # rateC = rate parameter of the exponential distribution of C
    rateC=1
    #rateC=15.5
    
    ######################################
    # parameters of inferential interest
    ######################################
    
    beta1 <- 0

    #beta2 <- matrix(c(0.3, -0.7, 0.1), nrow = p, byrow = T)
    #beta2 <- matrix(c(0.3, -0.7), nrow = p, byrow = T)
    beta2 <- matrix(0)
    
    #theta <- matrix(c(-0.6, -0.15, 0.4), nrow = p, byrow = T)
    #theta <- matrix(c(-0.6, -0.15), nrow = p, byrow = T)
    theta <- matrix(-0.2)
    
    ######################################
    # parameters not of inferential interest
    ######################################
    
    alpha <- matrix(rnorm(N*p, 0, 1), nrow = N, byrow = T)
    #alpha <- matrix(rgamma(N*p, 1, 1), nrow = N, byrow = T)
    #alpha <- matrix(log(rgamma(N*p, 1, 1)), nrow = N, byrow = T)
    #alpha <- matrix(as.numeric(scale(log(rgamma(N*p, 1, 1)))), nrow = N, byrow = T)
    #alpha <- matrix(as.numeric((log(rgamma(N*p, 1, 1))) - (digamma(1) - log(1)))/sqrt(trigamma(1)), nrow = N, byrow = T)
    #alpha <- matrix(runif(N*p, 0, sqrt(12)), nrow = N, byrow = T)
    #alpha <- matrix(rlogis(N*p, 0, sqrt(3)/pi), nrow = N, byrow = T)
    
    sigma <- 1
    
    ######################################
    # simulate data
    ######################################
    
    # m = number of occasions for all subject
    m <- rep(9, N)
    
    # max.m = max number of occasions
    max.m <- 10
    
    # A = treatment indicator
    A <- sample(c(0:1), N, replace=T, prob = c(0.3,0.7))
    #A <- sample(c(0:1), N, replace=T)
    
    # Generate survival outcome
    v <- runif(n=N)
    Tlat <- (- log(v) / (1* exp(t(t(theta) %*% t(alpha))+A*(t(beta1+t(beta2) %*% t(alpha))))))#^(1 / rho) # Weibull event times generated using IP methods
    #Tlat1 <- rexp(N, rate = exp(t(t(theta) %*% t(alpha))+A*(t(beta1+t(beta2) %*% t(alpha))))) # another way of generating exponential survival time
    C <- rexp(n=N, rate=rateC) # censoring times
    #end <- rep(0.105, N)
    end <- rep(0.251, N)
    time <- pmin(Tlat, C, end) # follow-up times
    #time <- pmin(Tlat, C)
    status <- as.numeric(Tlat == time) #event indicators
    dat <- data.frame(id=1:N,time=time,status=status, trt = A)
    
    # U = U_i for all subjects
    U <- numeric()
    for (j in 1:N){
      mi <- m[j]
      ui <- matrix(c(-9, -8, -7, -6, -5, -4, -3, -2, -1), nrow = mi)
      # ui <- matrix(c(1:mi), nrow = mi)
      Ui <- numeric()
      for (i in 1:p){
        Ui <- cbind(Ui, ui^(i-1))
      }
      extra.rows <- matrix(0, nrow = max.m-mi, ncol = p)
      Ui <- rbind(Ui, extra.rows)
      U <- rbind(U, Ui)
    }
    
    # X = longitudinal covariates
    alpha.comp <- alpha[rep(1:nrow(alpha), each = max.m), ]
    lin.pred <- matrix(rowSums(alpha.comp*U), nrow = N, byrow = T)
    e <- (lin.pred!=0)*matrix(rnorm(N*max.m, 0, sigma), nrow = N, byrow = T)
    X <- lin.pred+e
    dat <- cbind(dat, X)
    
    write.csv(dat,file=paste0(path, 'dat',z,'.csv'), row.names = F)
    write.csv(U,file=paste0(path, 'U', z,'.csv'), row.names = F)
    write.csv(alpha,file=paste0(path, 'alpha', z,'.csv'), row.names = F)
    
  })
}
  