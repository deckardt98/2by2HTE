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

#test1 <- data_gene(1,60,0.2,20,0.04,0.02,1,0.01074561)

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
    #t-test
    pvn[i] <- coef(summary(lmm1))[po,5]
    pva[i] <- coef(summary(lmm2))[po,5]
  }
  #calculate true power
  tp <- power_cal(type, nc, delta, m, rho0, rho1, sigma0, piz, pix)
  output <- as.data.frame(cbind(type, nc,round(mean(pvn[]<0.05, na.rm=T),3),round(tp,3),round(mean(pva[]<0.05, na.rm=T),3)))
  names(output) <- c("test.type","number.of.cluster","empirical.typeI.error","predicted.power","empirical.power")
  return(output)
}

B <- 5000
m <- c(20,50,100)
rho0 <- c(0.01,0.01,0.05,0.05,0.10,0.10)
rho1 <- c(0.02,0.04,0.07,0.10,0.12,0.15)
delta1 <- c(0.2,0.4)
delta2 <- c(0.15,0.3)

library(nlme)

#function generating simulation results
simre <- function(type, delta, m, rho0, rho1, sigma0, B = 5000){
  result <- c()
  for (j in 1:length(delta)){
    for (k in 1:length(m)){
      for (l in 1:length(rho0)){
        result <- rbind(result,simula(type, delta[j], m[k], rho0[l], rho1[l], sigma0, B = 5000))
      }
    }
  }
  return(result)
}

#function adding label
addla <- function(type,rho0,rho1,m,delta,resul){
  resul$type <- rep(type,36)
  resul$rho0 <- rep(rho0,6)
  resul$rho1 <- rep(rho1,6)
  resul$m <- rep(rep(m,c(6,6,6)),2)
  resul$delta <- rep(delta,c(18,18))
  return(resul)
}

rew1 <- simre(1,delta1,m,rho0,rho1,1,B=5000)
finalw1 <- addla(1,rho0,rho1,m,delta1,rew1)
write.csv(finalw1,"/Users/deckard/desktop/Fan Li Project/Project 2/t-test simulation results for EX without Interaction.csv")

rew2 <- simre(2,delta2,m,rho0,rho1,1,B=5000)
finalw2 <- addla(2,rho0,rho1,m,delta2,rew2)
write.csv(finalw2,"/Users/deckard/desktop/Fan Li Project/Project 2/t-test simulation results for EZ without Interaction.csv")
