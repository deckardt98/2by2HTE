#f(m)
f <- function(icc,m){
  return((icc*(m-1)+1))
}

#function calculating the number of cluster needed without interaction
nc <- function(type, delta, m, rho0, rho1, sigma0, piz, pix){
  #Input
  #type: test type 1=ex 2=ez
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #rho1: rho_1
  #sigma0: sigma_0^2
  #piz: pi_Z
  #pix: pi_X
  
  #output:
  #n: the number of cluster 
  
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #sigma_1^2
  sigma1 <- sigma0*(1-rho0)/(1-rho1)
    
  if (type==1){
    #variance
    v <- sigma0*f(rho0,m)/(1-pix)/m + sigma1*f(rho1,m)/pix/m
  }else if (type==2){
    #variance
    v <- 1/m/sz/((1-pix)*f(rho0,(m-1))/sigma0/(1-rho0)/f(rho0,m)+
                   pix*f(rho1,(m-1))/sigma1/(1-rho1)/f(rho1,m))
  }
  n <- (qnorm(1-alpha/2)+qnorm(1-mu))^2/delta^2*v
  return(n)
}

#Simulation Codes:

#Power calculation function without interaction
power_cal <- function(type, nc, delta, m, rho0, rho1, sigma0, piz, pix){
  
  #Input
  #type: test type: 1=ex 2=ez
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #rho1: rho_1
  #sigma0: sigma_0^2
  #piz: pi_Z
  #pix: pi_X
  
  
  #output:
  #power
  
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #sigma_1^2
  sigma1 <- sigma0*(1-rho0)/(1-rho1)
  
  if (type==1){
    #variance
    v <- sigma0*f(rho0,m)/(1-pix)/m + sigma1*f(rho1,m)/pix/m
  }else if (type==2){
    #variance
    v <- 1/m/sz/((1-pix)*f(rho0,(m-1))/sigma0/(1-rho0)/f(rho0,m)+
                   pix*f(rho1,(m-1))/sigma1/(1-rho1)/f(rho1,m))
  }

  power <- pnorm(sqrt(nc*delta^2/v) - qnorm(1-alpha/2))
  return(power)
}

#data-generation function
data_gene <- function(type, nc, delta, m, rho0, rho1, sigma0){
  
  #Input
  #type: test type: 1=ex 2=ez 
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #rho1: rho_1
  #sigma0: sigma_0^2
  
  
  #output:
  #sim_data: simulated data set
  
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  #sigma_1^2
  sigma1 <- sigma0*(1-rho0)/(1-rho1)
  
  #effect size
  beta1 = 1
  if (type==1){
    beta2 <- delta
    beta3 <- 0.05
  } else if (type==2){
    beta2 <- 0.15
    beta3 <- delta
  } 
  #cluster size
  cs <- rep(m, nc)
  #generate Xi
  ia <- rep(0, nc)
  rs <- sample(1:nc, size = nc*pix)
  ia[rs] = 1
  X <- rep(ia, cs)
  #generate Zij
  Z <- rbinom(sum(cs),1,piz)
  #generate outcome
  var_epsilon <- sigma0*(1-rho0)
  var_gamma <- sigma0-var_epsilon
  b2 <- sigma1-sigma0   #sigma_b_2^2
  gamma <- rep(rnorm(nc,0,sqrt(var_gamma)) ,cs)
  epsilon <- rnorm(sum(cs),0,sqrt(var_epsilon))
  sb2 <- rep(rnorm(nc,0,sqrt(b2)) ,cs)
  #alternative hypothesis
  Y_alt <- beta1 + (beta2+sb2)*X + beta3*Z + gamma + epsilon
  #null hypothesis
  if (type==1){
    Y_null <-  beta1 + sb2*X + beta3*Z + gamma + epsilon
  } else if (type==2){
    Y_null <- beta1 + (beta2+sb2)*X + gamma + epsilon
  } 
  #whole index
  ind <- rep(1:nc, cs)
  
  sim_data <- data.frame(ind, X, Z, Y_alt, Y_null)
  
  return(sim_data)
}

#function implementing simulation
simula <- function(type, delta, m, rho0, rho1, sigma0, B = 5000){
  
  #Input
  #type: test type: 1=ex 2=ez
  #delta: effect size
  #m: cluster size
  #rho0: rho_0
  #rho1: rho_1
  #sigma0: sigma_0^2
  #B: simulation times
  
  require(nlme)
  #type i error rate
  alpha <- 0.05
  #1-power
  mu <- 0.2
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  #sigma_1^2
  sigma1 <- sigma0*(1-rho0)/(1-rho1)
  
  #calculate the sample size via our formula
  mb <- m*piz
  nc <- ceiling(nc(type, delta, m, rho0, rho1, sigma0, piz, pix))
  #round up to the nearest even integer
  if (nc %% 2 ==1){nc <- nc + 1}
  #store p-values for null distribution
  #pvn <- array(NA,dim=B)
  pvn <- array(NA,dim=B)
  #store p-values for alternative distribution
  pva <- array(NA,dim=B)
  if (type==1){
    po <- 2
  } else if (type==2){
    po <- 3
  } 
  for (i in 1:B){
    
    #generate data
    #indicator for singular-fitting
    exit <- FALSE
    while (exit==FALSE){
      df <- data_gene(type, nc, delta, m, rho0, rho1, sigma0)
      lmm1 <- try(lme(Y_null ~  X+Z, data = df, 
                      random = list(ind = pdDiag(~ X))))
      lmm2 <- try(lme(Y_alt ~  X+Z, data = df, 
                      random = list(ind = pdDiag(~ X))))
      if(class(lmm1)!="try-error"&class(lmm2)!="try-error"){
        exit <- TRUE
      }
    }
    #z-test
    esn <- coef(summary(lmm1))[po,1]
    sdn <- coef(summary(lmm1))[po,2]
    testan <- esn/sdn
    pvn[i] <- (min((1-pnorm(testan)),pnorm(testan)))*2
    
    esa <- coef(summary(lmm2))[po,1]
    sda <- coef(summary(lmm2))[po,2]
    testaa <- esa/sda
    pva[i] <- (min((1-pnorm(testaa)),pnorm(testaa)))*2
  }
  #calculate true power
  tp <- power_cal(type, nc, delta, m, rho0, rho1, sigma0, piz, pix)
  output <- as.data.frame(cbind(type, 
                                delta,
                                m, 
                                rho0,
                                rho1,
                                nc,
                                round(mean(pvn[]<0.05, na.rm=T),3),
                                round(mean(pva[]<0.05, na.rm=T),3),
                                round(tp,3)))
  names(output) <- c("test.type",
                     "delta",
                     "m",
                     "alpha",
                     "alphaast",
                     "number.of.cluster",
                     "empirical.typeI.error",
                     "empirical.power",
                     "predicted.power")
  return(output)
}
library(parallel)
library(nlme)
B <- 5000
m <- c(20,50,100)
rho0 <- c(0.01,0.01,0.05,0.05,0.10,0.10)
rho1 <- c(0.02,0.04,0.07,0.10,0.12,0.15)
#delta1 <- c(0.15,0.3)
#delta2 <- c(0.1,0.15)
delta <- list(c(1,0.15),
              c(1,0.30),
              c(2,0.10),
              c(2,0.15))

#function generating simulation results
simre <- function(type, delta, m, rho0, rho1, sigma0, B = 5000, seed = 920784642){
  set.seed(seed)
  result <- c()
    for (k in 1:length(m)){
      for (l in 1:length(rho0)){
        result <- rbind(result,simula(type, delta, m[k], rho0[l], rho1[l], sigma0, B = B))
      }
    }
  return(result)
}
frt <- mclapply(delta, function(delta_val) {
  simre(delta_val[1],delta_val[2],m,rho0,rho1,1,B=5000, seed=920784642)
}, mc.cores = 4)
write.csv(frt[[1]],"/Users/deckard/desktop/Fan Li Project/Factorial HTE/z-test without interaction/new/ex0.15.csv")
write.csv(frt[[2]],"/Users/deckard/desktop/Fan Li Project/Factorial HTE/z-test without interaction/new/ex0.30.csv")
write.csv(frt[[3]],"/Users/deckard/desktop/Fan Li Project/Factorial HTE/z-test without interaction/new/ez0.10.csv")
write.csv(frt[[4]],"/Users/deckard/desktop/Fan Li Project/Factorial HTE/z-test without interaction/new/ez0.15.csv")
frt[[1]]
frt[[2]]
frt[[3]]
frt[[4]]
