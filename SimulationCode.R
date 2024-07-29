#########Simulation Code for the random slope model with intereaction###########

#Eigenvalue of the correlation matrix for the basic model, i.e., lambda(m) 
f <- function(icc,m){
  #Input:
  #icc: alpha or rho0
  #m: cluster size
  
  #Output: eigenvalue of the correlation matrix
  return((icc*(m-1)+1))
}

#function calculating series approximations for p0-p3
approxp <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_{b_4}^2
  #piz: pi_Z
  
  #Output:
  #app: a vector of (p0,p1,p2,p3)
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  #mean of m1i
  mb <- m*piz
  #sigma_Z^2
  sZ <- piz*(1-piz)
  
  #determinant
  D <- f(rho0,m)^2*b4^2+4*b4*rho0*ep*f(rho0,m)
  #roots
  x1 <- (f(rho0,m)*b4-sqrt(D))/(2*b4*rho0)
  x2 <- (f(rho0,m)*b4+sqrt(D))/(2*b4*rho0)
  #factors
  f1 <- ((x1+x2)*(x1+x2-3*mb)-x1*x2+3*mb^2)/((x1-mb)^3*(x2-mb)^3)
  f2 <- (x1+x2-2*mb)/((x1-mb)^2*(x2-mb)^2)-2*mb*f1
  f3 <- 3*mb^2*f1+(x1*x2-mb^2)/((x1-mb)^2*(x2-mb)^2)-
    2*mb*(x1*x2*(x1+x2-3*mb)+mb^3)/((x1-mb)^3*(x2-mb)^3)
  #third moment for m1i
  tm <- m*(m-1)*(m-2)*piz^3+3*m*(m-1)*piz^2+m*piz
  if (b4==0){
    app0 <- 1/(ep*f(rho0,m))
    app1 <- mb*app0
    app2 <- (m*sZ+mb^2)*app0
    app3 <- tm*app0
  }else {
    #approximated p0
    app0 <- -1/(b4*rho0)*(f1*(m*sZ+mb^2)+f2*mb+f3)
    #approximated p1
    app1 <- -1/(b4*rho0)*(f1*tm+f2*(m*sZ+mb^2)+f3*mb)
    #approximated p2
    app2 <- -1/(b4*rho0)+f(rho0,m)/rho0*app1+ep*f(rho0,m)/(b4*rho0)*app0
    #approximated p3
    app3 <- -mb/(b4*rho0)+f(rho0,m)/rho0*app2+ep*f(rho0,m)/(b4*rho0)*app1
  }
  #vector output
  app <- c(app0,app1,app2,app3)
  return(app)
}

#Function calculating lambda11 lambda12 lambda22 using series approximations
callamb <- function(m, rho0, alpha, sigma0, b4, piz){
  #Input
  #m: cluster size
  #rho0 
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  
  #output:
  #estlamb lambda11,lambda12,lambda22
  
  #sigma_epsilon^2
  ep <- sigma0*(1-alpha)
  esp <- approxp(m, rho0, alpha, sigma0, b4, piz)
  p0 <- esp[1]
  p1 <- esp[2]
  p2 <- esp[3]
  p3 <- esp[4]
  
  #calculate lambda11 lambda12 lambda22
  lamb11 <- (m+b4*(2*m*rho0-f(rho0,m))*p2-rho0*m^2*b4*p1-rho0*m^2*ep*p0)/(ep)
  lamb12 <- (m*piz+rho0*b4*p3-f(rho0,m)*b4*p2-m*rho0*ep*p1)/(ep)
  lamb22 <- (m*piz+rho0*b4*p3-(f(rho0,m)*b4+rho0*ep)*p2)/(ep)
  
  #output:
  estlamb <- c(lamb11,lamb12,lamb22)
  return(estlamb)
}

#Function calculating the number of cluster using series approximations
nc <- function(type, delta, m, rho0, alpha, sigma0, b4, piz, pix){
  #Input
  #type: test type 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  #Output:
  #n: the number of cluster 
  
  #type I error rate
  epsilon1 <- 0.05
  #1-power
  epsilon2 <- 0.2
  #sigma_Z^2
  sz <- piz*(1-piz)
  #calculate lambda
  lamb <- callamb(m, rho0, alpha, sigma0, b4, piz)
  l11 <- lamb[1]
  l12 <- lamb[2]
  l22 <- lamb[3]
  if (type==1){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*((m-1)*alpha*(1-piz)+1-alpha))/(m*(1-pix)*(1-piz)*((m-2)*alpha+1))
    #variance part2
    vp2 <- l22/(pix*(l11*l22-l12^2))
  }else if (type==2){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*(1-alpha))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- 0
  }else if (type==3){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1))/(m*(1-pix))
    #variance part2
    vp2 <- (piz^2*l11-2*piz*l12+l22)/(pix*(l11*l22-l12^2))
  } else if (type==4){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1)*(1-pix))/(m*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- (pix*l11)/(l11*l22-l12^2)
  }else if (type==5){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- l11/(pix*(l11*l22-l12^2))
  }
  n <- (qnorm(1-epsilon1/2)+qnorm(1-epsilon2))^2/delta^2*(vp1+vp2)
  return(n)
}

#Function for power calculation  
power_cal <- function(type, nc, delta, m, rho0, alpha, sigma0, b4, piz, pix){
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #piz: pi_Z
  #pix: pi_X
  
  
  #Output:
  #power
  
  #type i error rate
  epsilon1 <- 0.05
  #sigma_Z^2
  sz <- piz*(1-piz)
  #calculate lambda
  lamb <- callamb(m, rho0, alpha, sigma0, b4, piz)
  l11 <- lamb[1]
  l12 <- lamb[2]
  l22 <- lamb[3]
  if (type==1){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*((m-1)*alpha*(1-piz)+1-alpha))/(m*(1-pix)*(1-piz)*((m-2)*alpha+1))
    #variance part2
    vp2 <- l22/(pix*(l11*l22-l12^2))
  }else if (type==2){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1)*(1-alpha))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- 0
  }else if (type==3){
    #variance part 1
    vp1 <- (sigma0*((m-1)*alpha+1))/(m*(1-pix))
    #variance part2
    vp2 <- (piz^2*l11-2*piz*l12+l22)/(pix*(l11*l22-l12^2))
  } else if (type==4){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1)*(1-pix))/(m*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- (pix*l11)/(l11*l22-l12^2)
  }else if (type==5){
    #variance part 1
    vp1 <- (sigma0*(1-alpha)*((m-1)*alpha+1))/(m*(1-pix)*sz*((m-2)*alpha+1))
    #variance part2
    vp2 <- l11/(pix*(l11*l22-l12^2))
  }
  var <- vp1+vp2
  power <- pnorm(sqrt(nc*delta^2/var) - qnorm(1-epsilon1/2))
  return(power)
}

#data-generation function
data_gene <- function(type, nc, delta, m, rho0, alpha, sigma0, b4){
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #nc: number of clusters
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  
  
  #output:
  #sim_data: simulated data set
  
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  
  #effect size
  beta1 = 1
  if (type==1){
    beta2 <- delta
    beta3 <- 0.05
    beta4 <- 0.05
  } else if (type==2){
    beta2 <- 0.15
    beta3 <- delta
    beta4 <- 0.05
  } else if (type==3){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 2*delta-0.3
  } else if (type==4){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- 2*delta-0.1
  } else if (type==5){
    beta2 <- 0.15
    beta3 <- 0.05
    beta4 <- delta
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
  var_epsilon <- sigma0*(1-alpha)
  var_gamma <- sigma0-var_epsilon
  b2 <- sigma0*(rho0-alpha)/(1-rho0)   #sigma_b_2^2
  gamma <- rep(rnorm(nc,0,sqrt(var_gamma)) ,cs)
  epsilon <- rnorm(sum(cs),0,sqrt(var_epsilon))
  sb2 <- rep(rnorm(nc,0,sqrt(b2)) ,cs)
  sb4 <- rep(rnorm(nc,0,sqrt(b4)) ,cs)
  #alternative hypothesis
  Y_alt <- beta1 + (beta2+sb2)*X + beta3*Z + (beta4+sb4)*X*Z + gamma + epsilon
  #null hypothesis
  if (type==1){
    Y_null <-  beta1 + sb2*X + beta3*Z + (beta4+sb4)*X*Z + gamma + epsilon
  } else if (type==2){
    Y_null <- beta1 + (beta2+sb2)*X + (beta4+sb4)*X*Z + gamma + epsilon
  } else if (type==3){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + (-0.3+sb4)*X*Z + gamma + epsilon
  } else if (type==4){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + (-0.1+sb4)*X*Z + gamma + epsilon
  } else if (type==5){
    Y_null <- beta1 + (beta2+sb2)*X + beta3*Z + sb4*X*Z + gamma + epsilon
  }
  
  #whole index
  ind <- rep(1:nc, cs)
  
  sim_data <- data.frame(ind, X, Z, Y_alt, Y_null)
  
  return(sim_data)
}

#Function that implements one single scenario of the simulation
simula <- function(type, delta, m, rho0, alpha, sigma0, b4, rho1, B = 5000){
  #Input
  #type: test type: 1=cex 2=cez 3=mex 4=mez 5=ie
  #delta: effect size
  #m: cluster size
  #rho0
  #alpha
  #sigma0: sigma_0^2
  #b4: sigma_b_4^2
  #rho1
  #B: the number of iterations per scenario  
  
  require(nlme)
  #type i error rate
  epsilon1 <- 0.05
  #pi_X
  pix <- 0.5
  #pi_Z
  piz <- 0.5
  
  #calculate the sample size via our formula
  mb <- m*piz
  nc <- ceiling(nc(type, delta, m, rho0, alpha, sigma0, b4, piz, pix))
  #round up to the nearest even integer
  if (nc %% 2 ==1){nc <- nc + 1}
  #store p-values for null distribution
  #pvn <- array(NA,dim=B)
  pvn <- array(NA,dim=B)
  #store p-values for alternative distribution
  pva <- array(NA,dim=B)
  #position
  po2 <- 4
  if (type==1){
    po1 <- 2
    co <- 0
  } else if (type==2){
    po1 <- 3
    co <- 0
  } else if (type==3){
    po1 <- 2
    co <- piz
  } else if (type==4){
    po1 <- 3
    co <- pix
  } else if (type==5){
    po1 <- 4
    co <- 0
  }
  for (i in 1:B){
    
    #generate data
    #indicator for singular-fitting
    exit <- FALSE
    while (exit==FALSE){
      df <- data_gene(type, nc, delta, m, rho0, alpha, sigma0, b4)
      lmm1 <- try(lme(Y_null ~  X*Z, data = df, 
                      random = list(ind = pdDiag(~ X*Z-Z))))
      lmm2 <- try(lme(Y_alt ~  X*Z, data = df, 
                      random = list(ind = pdDiag(~ X*Z-Z))))
      if(class(lmm1)!="try-error"&class(lmm2)!="try-error"){
        exit <- TRUE
      }
    }
    #z-test
    esn <- coef(summary(lmm1))[po1,1]+co*coef(summary(lmm1))[po2,1]
    esvn <- coef(summary(lmm1))[po1,2]^2+co^2*coef(summary(lmm1))[po2,2]^2+2*co*vcov(lmm1)[po1,po2]
    testan <- esn/sqrt(esvn)
    pvn[i] <- (min((1-pnorm(testan)),pnorm(testan)))*2
    
    esa <- coef(summary(lmm2))[po1,1]+co*coef(summary(lmm2))[po2,1]
    esva <- coef(summary(lmm2))[po1,2]^2+co^2*coef(summary(lmm2))[po2,2]^2+2*co*vcov(lmm2)[po1,po2]
    testaa <- esa/sqrt(esva)
    pva[i] <- (min((1-pnorm(testaa)),pnorm(testaa)))*2
  }
  #calculate true power
  tp <- power_cal(type, nc, delta, m, rho0, alpha, sigma0, b4, piz, pix)
  output <- as.data.frame(cbind(type, 
                                delta,
                                m, 
                                alpha,
                                rho0,
                                rho1,
                                nc,round(mean(pvn[]<0.05, na.rm=T),3),
                                round(mean(pva[]<0.05, na.rm=T),3),
                                round(tp,3)))
  names(output) <- c("test.type",
                     "delta",
                     "m",
                     "alpha",
                     "rho0",
                     "rho1",
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
alpha <- c(0.01,0.01,0.01,0.01,0.05,0.05,0.10,0.10)
rho0 <- c(0.02,0.02,0.04,0.04,0.06,0.10,0.12,0.13)
rho1 <- c(0.05,0.12,0.05,0.12,0.08,0.13,0.14,0.15)

alpha <-  c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.10,0.10,0.10)
rho0 <- c(0.01,0.01,0.01,0.01,0.05,0.05,0.05,0.10,0.10,0.20,0.10,0.10,0.20)
rho1 <- c(0.01,0.05,0.10,0.20,0.05,0.10,0.20,0.10,0.20,0.20,0.10,0.20,0.20)
#alpha <- c(0.05,0.05)
#rho0 <- c(0.06,0.10)
#rho1 <- c(0.08,0.13)
b4 <- c()
for (i in 1:length(alpha)){
  b4[i] <- (1-alpha[i])*(rho1[i]-rho0[i])/(1-rho0[i])/(1-rho1[i])
}
delta <- list(c(1,0.2),
              c(1,0.35),
              c(2,0.15),
              c(2,0.3),
              c(3,0.2),
              c(3,0.35),
              c(4,0.1),
              c(4,0.15),
              c(5,0.2),
              c(5,0.3))
delta <- list(c(1,0.30),
              c(2,0.25),
              c(3,0.30),
              c(4,0.18),
              c(5,0.30))
#delta1 <- c(0.2,0.35)
#delta2 <- c(0.15,0.3)
#delta3 <- c(0.2,0.35)
#delta4 <- c(0.1,0.15)
#delta5 <- c(0.2,0.3)
#function generating simulation results
simre <- function(type, delta, m, rho0, alpha, sigma0, b4, rho1, B = 5000, seed = 920784642){
  set.seed(seed)
  result <- c()
  for (k in 1:length(m)){
    for (l in 1:length(rho0)){
      result <- rbind(result,simula(type, delta, m[k], rho0[l], alpha[l], sigma0, b4[l], rho1[l], B = B))
      #result <- c(result,nc(type, delta, m[k], rho0[l], alpha[l], sigma0, b4[l], 0.5, 0.5))
    }
  }
  return(result)
}
frt <- mclapply(delta, function(delta_val) {
  simre(delta_val[1],delta_val[2],m,rho0,alpha,1,b4,rho1,B=5, seed=920784642)
}, mc.cores = 5)
table3 <- cbind(frt[[2]][,-c(1,2)],frt[[4]][,-c(1:6)])
table4 <- cbind(frt[[6]][,-c(1,2)],frt[[8]][,-c(1:6)],frt[[10]][,-c(1:6)])
webtable3 <- cbind(frt[[1]][,-c(1,2)],frt[[3]][,-c(1:6)])
webtable4 <- cbind(frt[[5]][,-c(1,2)],frt[[7]][,-c(1:6)],frt[[9]][,-c(1:6)])
setwd("/Users/deckard/desktop/Fan Li Project/Factorial HTE/z-test with interaction/new")
write.csv(table3, "./table3.csv")
write.csv(table4, "./table4.csv")
write.csv(webtable3, "./webtable3.csv")
write.csv(webtable4, "./webtable4.csv")
